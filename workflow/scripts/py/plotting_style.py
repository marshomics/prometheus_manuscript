"""Shared plotting style for publication-ready, Adobe-Illustrator-editable figures.

Every Python figure imports `apply_style()` and saves through `save_fig()`, which
writes BOTH a PNG (high-DPI raster) and an SVG (vector). Font handling is chosen so
the output opens in Illustrator with TEXT THAT STAYS EDITABLE:

  * svg.fonttype = 'none'  -> text is written as <text> elements referencing the
    font family (editable as text in Illustrator; not converted to paths).
  * pdf.fonttype / ps.fonttype = 42 -> if a PDF/EPS is ever exported, TrueType
    fonts are embedded (Type 42), again keeping text editable.
  * A standard sans family (Arial, falling back to DejaVu Sans which ships with
    matplotlib) is used so glyphs are present on any machine.

Note on "embedded" SVG fonts: matplotlib does not base64-embed font files in SVG.
The defensible, standard choice for an Illustrator workflow is editable <text>
referencing a common family (done here). If a teammate lacks Arial, Illustrator
substitutes the family without rasterizing the text. For guaranteed-identical
rendering on any machine, also keep the PNG (raster) which is always produced.
"""
from __future__ import annotations
import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

DPI = int(os.environ.get("PROM_DPI", "600"))
FONT_FAMILY = os.environ.get("PROM_FONT", "Arial")
FALLBACK_FONT = os.environ.get("PROM_FALLBACK_FONT", "DejaVu Sans")
BASE_SIZE = float(os.environ.get("PROM_FONT_SIZE", "7"))

# Colour-blind-safe qualitative palette (Okabe–Ito).
OKABE_ITO = {
    "black": "#000000", "orange": "#E69F00", "skyblue": "#56B4E9",
    "green": "#009E73", "yellow": "#F0E442", "blue": "#0072B2",
    "vermillion": "#D55E00", "purple": "#CC79A7", "grey": "#999999",
}
TIER_COLORS = {"A": "#0072B2", "B": "#56B4E9", "C": "#999999", "Rejected": "#D55E00"}


def _resolve_family() -> str:
    available = {f.name for f in font_manager.fontManager.ttflist}
    return FONT_FAMILY if FONT_FAMILY in available else FALLBACK_FONT


def apply_style(font_family: str | None = None, base_size: float | None = None) -> None:
    fam = font_family or _resolve_family()
    size = base_size if base_size is not None else BASE_SIZE
    plt.rcParams.update({
        # ---- Illustrator-editable text ----
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        # ---- fonts ----
        "font.family": "sans-serif",
        "font.sans-serif": [fam, FALLBACK_FONT, "Helvetica", "Arial"],
        "font.size": size,
        "axes.titlesize": size + 1,
        "axes.labelsize": size,
        "xtick.labelsize": size - 1,
        "ytick.labelsize": size - 1,
        "legend.fontsize": size - 1,
        "figure.titlesize": size + 2,
        # ---- clean, journal-style axes ----
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "lines.linewidth": 1.0,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "figure.dpi": 150,
        "pdf.compression": 6,
    })


def default_formats() -> tuple:
    """config plotting.formats as exported by the workflow (PROM_PLOT_FORMATS, set in
    common.smk), else png + svg. The figure rules declare one output per listed format."""
    env = os.environ.get("PROM_PLOT_FORMATS", "")
    fm = tuple(x.strip().lower().lstrip(".") for x in env.split(",") if x.strip())
    return fm or ("png", "svg")


def _formats(formats) -> tuple:
    if formats is None:
        return default_formats()
    if isinstance(formats, str):
        formats = formats.split(",")
    fm = tuple(str(x).strip().lower().lstrip(".") for x in formats if str(x).strip())
    return fm or default_formats()


def save_fig(fig, prefix: str, formats=None, dpi: int | None = None) -> None:
    """Write `prefix.<ext>` for every format (default: config plotting.formats via
    PROM_PLOT_FORMATS, else png + svg). Raises if a file is not written."""
    os.makedirs(os.path.dirname(os.path.abspath(prefix)), exist_ok=True)
    d = dpi or DPI
    try:
        for ext in _formats(formats):
            path = f"{prefix}.{ext}"
            fig.savefig(path, dpi=d, transparent=True)
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                raise RuntimeError(f"figure not written: {path}")
    finally:
        plt.close(fig)


def placeholder(prefix: str, message: str, formats=None) -> None:
    """Emit a labelled empty figure so the DAG completes when an input is empty."""
    apply_style()
    fig, ax = plt.subplots(figsize=(4, 2.2))
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True, fontsize=BASE_SIZE)
    ax.axis("off")
    save_fig(fig, prefix, formats)
