# =====================================================================
# Shared plotting helpers for the R figure scripts.
#
# Every R figure script sources this file and saves through save_plot(),
# which writes BOTH:
#   * {prefix}.svg via svglite::svglite -- text is written as EDITABLE
#     <text> elements (not converted to paths/outlines), so the figure
#     opens in Adobe Illustrator with type that can still be selected,
#     re-flowed and re-styled. svglite also embeds the font glyphs needed
#     for faithful on-screen rendering while keeping the text editable.
#   * {prefix}.png via ragg::agg_png at res=600 -- a high-resolution
#     raster that renders identically on any machine regardless of which
#     fonts are installed (the safe fallback for figure review / journals
#     that want a flat raster).
#
# Font handling: we prefer Arial through the systemfonts registry so the
# label family matches the Python matplotlib figures (plotting_style.py).
# If Arial is not registered on the host, we fall back to the generic
# "sans" family, which every backend resolves. No call ever errors merely
# because a font is missing.
#
# This file defines only functions + constants; sourcing it has no side
# effects beyond (best-effort, guarded) loading of ggplot2.
# =====================================================================

suppressWarnings(suppressMessages({
  ok_gg <- requireNamespace("ggplot2", quietly = TRUE)
}))

# ---- colour-blind-safe qualitative palette (Okabe-Ito) ---------------
# Matches plotting_style.py so Python and R panels share one palette.
OKABE_ITO <- c(
  black      = "#000000", orange   = "#E69F00", skyblue = "#56B4E9",
  green      = "#009E73", yellow   = "#F0E442", blue    = "#0072B2",
  vermillion = "#D55E00", purple   = "#CC79A7", grey    = "#999999"
)

# Confidence-tier colours (A/B/C/Rejected) shared with the Python figures.
TIER_COLORS <- c(A = "#0072B2", B = "#56B4E9", C = "#999999", Rejected = "#D55E00")

# Presence / methanogen two-level fills.
PRESENCE_COLORS <- c(`TRUE` = "#0072B2", `FALSE` = "#D9D9D9",
                     present = "#0072B2", absent = "#D9D9D9")
METHANOGEN_COLORS <- c(`TRUE` = "#009E73", `FALSE` = "#F0E442")


# --------------------------------------------------------------------- #
# Pick a usable font family. Prefer Arial via systemfonts; else "sans".
# Returns a length-1 character; never errors.
# --------------------------------------------------------------------- #
.pick_font_family <- function(preferred = "Arial") {
  fam <- "sans"
  tryCatch({
    if (requireNamespace("systemfonts", quietly = TRUE)) {
      faces <- systemfonts::system_fonts()
      if (is.data.frame(faces) && "family" %in% names(faces) &&
          preferred %in% faces$family) {
        fam <- preferred
      }
    }
  }, error = function(e) invisible(NULL))
  fam
}

PUB_FONT <- .pick_font_family("Arial")


# --------------------------------------------------------------------- #
# theme_pub(base_size = 7)
# A clean, journal-style ggplot2 theme: white background, thin axes,
# no top/right spines, small legend, Arial (or sans fallback) family.
# --------------------------------------------------------------------- #
theme_pub <- function(base_size = 7, base_family = PUB_FONT) {
  if (!isTRUE(ok_gg)) return(NULL)  # ggplot2 missing -> no-op (guarded use)
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "grey90"),
      panel.border     = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(linewidth = 0.4, colour = "black"),
      axis.ticks       = ggplot2::element_line(linewidth = 0.4, colour = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(size = base_size, face = "bold"),
      legend.key       = ggplot2::element_blank(),
      legend.key.size  = ggplot2::unit(0.9, "lines"),
      legend.title     = ggplot2::element_text(size = base_size),
      legend.text      = ggplot2::element_text(size = base_size - 1),
      plot.title       = ggplot2::element_text(size = base_size + 1, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = base_size - 1, colour = "grey30"),
      complete = TRUE
    )
}


# --------------------------------------------------------------------- #
# save_plot(plot, prefix, width = 4, height = 3)
# Write {prefix}.png (ragg::agg_png, res=600 raster) AND {prefix}.svg
# (svglite::svglite, EDITABLE <text> for Illustrator). `plot` may be a
# ggplot, a patchwork/cowplot object, or a function/expression that draws
# to the active device (for base-graphics figures such as cophyloplot /
# phytools maps). Dimensions are in inches. Never throws: on any device
# error it writes a minimal placeholder so the DAG still completes.
# --------------------------------------------------------------------- #
save_plot <- function(plot, prefix, width = 4, height = 3, res = 600) {
  dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
  is_grob <- isTRUE(ok_gg) && (ggplot2::is.ggplot(plot) ||
                               inherits(plot, c("patchwork", "gg", "ggassemble")))
  is_drawfun <- is.function(plot)

  .render <- function(dev_open, dev_close) {
    opened <- FALSE
    tryCatch({
      dev_open(); opened <- TRUE
      if (is_grob) {
        print(plot)
      } else if (is_drawfun) {
        plot()
      } else if (inherits(plot, "recordedplot")) {
        grDevices::replayPlot(plot)
      } else {
        # last resort: try printing whatever it is
        print(plot)
      }
    }, error = function(e) {
      message("save_plot: render error: ", conditionMessage(e))
      tryCatch({
        graphics::plot.new()
        graphics::text(0.5, 0.5, paste("Figure unavailable:", conditionMessage(e)))
      }, error = function(e2) invisible(NULL))
    }, finally = {
      if (opened) tryCatch(dev_close(), error = function(e) invisible(NULL))
    })
  }

  png_path <- paste0(prefix, ".png")
  svg_path <- paste0(prefix, ".svg")

  # ---- PNG (ragg high-res raster) ---- #
  if (requireNamespace("ragg", quietly = TRUE)) {
    .render(function() ragg::agg_png(png_path, width = width, height = height,
                                     units = "in", res = res, background = "white"),
            grDevices::dev.off)
  } else {
    # Fallback raster if ragg unavailable.
    .render(function() grDevices::png(png_path, width = width, height = height,
                                      units = "in", res = res, bg = "white"),
            grDevices::dev.off)
  }

  # ---- SVG (svglite, editable text) ---- #
  if (requireNamespace("svglite", quietly = TRUE)) {
    .render(function() svglite::svglite(svg_path, width = width, height = height,
                                        bg = "white"),
            grDevices::dev.off)
  } else {
    .render(function() grDevices::svg(svg_path, width = width, height = height,
                                      bg = "white"),
            grDevices::dev.off)
  }

  invisible(c(png = png_path, svg = svg_path))
}


# --------------------------------------------------------------------- #
# placeholder_plot(prefix, message)
# A labelled empty figure (png + svg) for trivial/empty inputs, so the
# Snakemake DAG always finds its declared outputs. Uses save_plot.
# --------------------------------------------------------------------- #
placeholder_plot <- function(prefix, message, width = 4, height = 2.2) {
  draw <- function() {
    op <- graphics::par(mar = c(0, 0, 0, 0))
    on.exit(graphics::par(op), add = TRUE)
    graphics::plot.new()
    graphics::text(0.5, 0.5, message, cex = 0.9, col = "grey20")
  }
  save_plot(draw, prefix, width = width, height = height)
}
