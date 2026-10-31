# =====================================================================
# Shared plotting helpers for the R figure scripts.
#
# Every R figure script sources this file and saves through save_plot(),
# which writes every format in config plotting.formats (default png + svg):
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

PUB_FONT <- .pick_font_family(Sys.getenv("PROM_FONT", "Arial"))
# config plotting.* as exported by common.smk (PROM_*); defaults = the documented ones.
PUB_BASE_SIZE <- suppressWarnings(as.numeric(Sys.getenv("PROM_FONT_SIZE", "7")))
if (is.na(PUB_BASE_SIZE) || PUB_BASE_SIZE <= 0) PUB_BASE_SIZE <- 7
PLOT_DPI <- suppressWarnings(as.integer(Sys.getenv("PROM_DPI", "600")))
if (is.na(PLOT_DPI) || PLOT_DPI < 72) PLOT_DPI <- 600L

# Formats every figure is written in: config plotting.formats (PROM_PLOT_FORMATS), the
# same list the figure rules declare as outputs. Supported: png, svg, pdf.
plot_formats <- function() {
  f <- tolower(sub("^\\.", "", trimws(strsplit(Sys.getenv("PROM_PLOT_FORMATS", "png,svg"), ",")[[1]])))
  f <- f[nzchar(f)]
  if (!length(f)) f <- c("png", "svg")
  bad <- setdiff(f, c("png", "svg", "pdf"))
  if (length(bad)) stop("plotting.formats: unsupported format(s): ", paste(bad, collapse = ", "))
  unique(f)
}


# --------------------------------------------------------------------- #
# theme_pub(base_size = 7)
# A clean, journal-style ggplot2 theme: white background, thin axes,
# no top/right spines, small legend, Arial (or sans fallback) family.
# --------------------------------------------------------------------- #
theme_pub <- function(base_size = PUB_BASE_SIZE, base_family = PUB_FONT) {
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
# Write {prefix}.<ext> for every configured format (plot_formats(): png via
# ragg::agg_png at PLOT_DPI, svg via svglite with EDITABLE <text> for Illustrator,
# pdf via cairo_pdf). `plot` may be a ggplot, a patchwork/cowplot object, or a
# function that draws to the active device (base graphics: cophyloplot, phytools
# maps). Dimensions are in inches. A render error FAILS (the partial file is removed
# and the error propagates), so a broken figure can never pass as a finished output;
# placeholders are drawn only deliberately, via placeholder_plot() for empty inputs.
# --------------------------------------------------------------------- #
save_plot <- function(plot, prefix, width = 4, height = 3, res = PLOT_DPI,
                      formats = plot_formats()) {
  dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
  is_grob <- isTRUE(ok_gg) && (ggplot2::is.ggplot(plot) ||
                               inherits(plot, c("patchwork", "gg", "ggassemble")))
  is_drawfun <- is.function(plot)

  .render <- function(path, dev_open) {
    opened <- FALSE
    ok <- FALSE
    err <- ""
    tryCatch({
      dev_open(); opened <- TRUE
      if (is_grob) {
        print(plot)
      } else if (is_drawfun) {
        plot()
      } else if (inherits(plot, "recordedplot")) {
        grDevices::replayPlot(plot)
      } else {
        print(plot)
      }
      ok <- TRUE
    }, error = function(e) {
      err <<- conditionMessage(e)
      message("save_plot: render error: ", err)
    }, finally = {
      if (opened) tryCatch(grDevices::dev.off(), error = function(e) invisible(NULL))
    })
    if (!ok) {
      unlink(path)
      stop("save_plot: failed to render ", path, ": ", err, call. = FALSE)
    }
    if (!file.exists(path) || isTRUE(file.info(path)$size == 0))
      stop("save_plot: figure not written: ", path, call. = FALSE)
    invisible(path)
  }

  paths <- character(0)
  for (fmt in formats) {
    path <- paste0(prefix, ".", fmt)
    opener <- switch(fmt,
      png = if (requireNamespace("ragg", quietly = TRUE)) {
        function() ragg::agg_png(path, width = width, height = height, units = "in",
                                 res = res, background = "white")
      } else {
        function() grDevices::png(path, width = width, height = height, units = "in",
                                  res = res, bg = "white")
      },
      svg = if (requireNamespace("svglite", quietly = TRUE)) {
        function() svglite::svglite(path, width = width, height = height, bg = "white")
      } else {
        function() grDevices::svg(path, width = width, height = height, bg = "white")
      },
      pdf = if (isTRUE(capabilities("cairo"))) {
        function() grDevices::cairo_pdf(path, width = width, height = height, bg = "white")
      } else {
        function() grDevices::pdf(path, width = width, height = height, bg = "white",
                                  useDingbats = FALSE)
      },
      stop("save_plot: unsupported format '", fmt, "'"))
    .render(path, opener)
    paths[[fmt]] <- path
  }
  invisible(paths)
}


# --------------------------------------------------------------------- #
# placeholder_plot(prefix, message)
# A labelled empty figure (configured formats) for trivial/empty inputs, so the
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
