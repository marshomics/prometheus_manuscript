"""Figure for the cross-domain Prometheus tree: tips coloured by domain + summary.

Draws the combined tree (Bio.Phylo) with tips coloured Archaea / Bacteria / Reference
when the tree is small enough to label, and always prints the headline metrics
(independent archaeal acquisitions, cross-domain transitions, inferred direction) in
the title. Writes PNG (600 dpi) + SVG with editable text. Any failure -> a placeholder
panel so the DAG never breaks. Pure plotting; no statistics here.
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

COLORS = {"Archaea": "#0072B2", "Bacteria": "#D55E00", "Reference": "#999999",
          "Paralog": "#CCCCCC", "Unknown": "#000000"}


def load_tips(path):
    dom = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    dom[f[0]] = f[1]
    return dom


def load_summary(path):
    s = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    s[f[0]] = f[1]
    return s


def domain_class(d):
    if d in ("Archaea", "Bacteria"):
        return d
    if d in ("Paralog",):
        return "Paralog"
    return "Reference" if d != "Unknown" else "Unknown"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--tips", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)), exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams["svg.fonttype"] = "none"
    import matplotlib.pyplot as plt

    dom = load_tips(args.tips)
    summ = load_summary(args.summary)
    title = "Cross-domain Prometheus"
    if summ:
        title += (f"  |  archaeal acquisitions={summ.get('n_archaeal_clades','?')}, "
                  f"transitions={summ.get('n_crossdomain_transitions','?')}, "
                  f"{summ.get('inferred_direction','?')}")

    def placeholder(msg):
        for ext in ("png", "svg"):
            fig, ax = plt.subplots(figsize=(6, 3))
            ax.axis("off")
            ax.text(0.5, 0.5, msg, ha="center", va="center", wrap=True, fontsize=8)
            fig.savefig(f"{args.out_prefix}.{ext}", dpi=600 if ext == "png" else None,
                        bbox_inches="tight")
            plt.close(fig)

    txt = ""
    if os.path.exists(args.tree) and os.path.getsize(args.tree) > 0:
        with open(args.tree) as fh:
            txt = fh.read().strip()
    if not txt or txt in ("();", ";") or txt.count(",") < 2:
        placeholder(title + "\n(no tree)")
        utils.log("plot_crossdomain_tree: trivial tree; placeholder.")
        return

    try:
        from io import StringIO
        from Bio import Phylo
        tree = Phylo.read(StringIO(txt), "newick")
        n_tips = tree.count_terminals()
        show_labels = n_tips <= 120

        def label_color(clade):
            name = getattr(clade, "name", None)
            return COLORS.get(dom.get(name, "Unknown"), "#000000")

        height = max(3.0, min(0.16 * n_tips, 24))
        fig, ax = plt.subplots(figsize=(8, height))
        Phylo.draw(
            tree, axes=ax, do_show=False,
            label_func=(lambda c: (c.name if (c.is_terminal() and show_labels) else "")),
            label_colors=label_color,
        )
        ax.set_title(title, fontsize=8)
        # legend
        from matplotlib.lines import Line2D
        handles = [Line2D([0], [0], marker="o", color="w", label=k,
                          markerfacecolor=v, markersize=6)
                   for k, v in (("Archaea", COLORS["Archaea"]),
                                ("Bacteria", COLORS["Bacteria"]),
                                ("Reference", COLORS["Reference"]))]
        ax.legend(handles=handles, fontsize=7, loc="lower right", frameon=False)
        for ext in ("png", "svg"):
            fig.savefig(f"{args.out_prefix}.{ext}", dpi=600 if ext == "png" else None,
                        bbox_inches="tight")
        plt.close(fig)
        utils.log(f"plot_crossdomain_tree: drew {n_tips} tips (labels={show_labels}).")
    except Exception as exc:
        placeholder(title + f"\n(render failed: {exc})")
        utils.log(f"plot_crossdomain_tree: render failed ({exc}); placeholder.")


if __name__ == "__main__":
    main()
