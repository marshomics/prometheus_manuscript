"""Figure for the cross-domain Prometheus tree: branches coloured by domain of life.

Draws the OUTGROUP-ROOTED tree written by crossdomain_analysis.py (--out-rooted-tree)
as a rectangular phylogram. Each branch is coloured by the composition of the tips
below it (all Archaea / all Bacteria / paralog outgroup / mixed), tips carry a
domain-coloured marker (so large trees stay readable without labels), nodes whose
SH-aLRT/UFBoot reach the support thresholds get a black dot, and the headline metrics
(robust and ML-tree direction, archaeal clades, transitions) go in the title. Tip
labels are drawn only for <= --max-labels tips. Writes every --formats extension (SVG
text stays editable). A trivial tree (no cross-domain analysis possible) gives a
labelled placeholder; any other failure is an error. Pure plotting; no statistics.
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO  # noqa: E402
from crossdomain_analysis import parse_newick, parse_support  # noqa: E402
import utils  # noqa: E402

COLORS = {"Archaea": OKABE_ITO["blue"], "Bacteria": OKABE_ITO["vermillion"],
          "Reference": OKABE_ITO["orange"], "Paralog": OKABE_ITO["grey"],
          "mixed": "#4D4D4D"}


def load_tips(path):
    dom, src = {}, {}
    with open(path) as fh:
        header = next(fh, "").rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(header, line.rstrip("\n").split("\t")))
            if f.get("tip"):
                dom[f["tip"]] = f.get("domain", "")
                src[f["tip"]] = f.get("source", "")
    return dom, src


def load_summary(path):
    s = {}
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                s[f[0]] = f[1]
    return s


def tip_class(lbl, dom, src):
    if src.get(lbl) == "prometheus_reference":
        return "Reference"
    d = dom.get(lbl, "")
    return d if d in ("Archaea", "Bacteria", "Paralog") else "mixed"


def layout(parent, children, length):
    """x = root distance (depth when no lengths), y = leaf order; returns x, y, order."""
    n = len(parent)
    use_len = any(v is not None and v > 0 for v in length[1:])
    x, order, stack = [0.0] * n, [], [0]
    while stack:
        v = stack.pop()
        order.append(v)
        for c in reversed(children[v]):
            x[c] = x[v] + ((length[c] or 0.0) if use_len else 1.0)
            stack.append(c)
    y, k = [0.0] * n, 0
    for v in order:
        if not children[v]:
            y[v] = k
            k += 1
    for v in reversed(order):
        if children[v]:
            ys = [y[c] for c in children[v]]
            y[v] = (min(ys) + max(ys)) / 2.0
    return x, y, order


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True, help="outgroup-rooted tree (crossdomain_analysis)")
    ap.add_argument("--tips", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--formats", default="png,svg")
    ap.add_argument("--max-labels", type=int, default=120)
    ap.add_argument("--min-ufboot", type=float, default=95.0)
    ap.add_argument("--min-alrt", type=float, default=80.0)
    args = ap.parse_args()
    formats = [f.strip() for f in args.formats.split(",") if f.strip()]

    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)), exist_ok=True)
    apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    summ = load_summary(args.summary)
    title = "Cross-domain Prometheus"
    if summ.get("status") == "complete":
        title += (f" | direction: {summ.get('inferred_direction', '?')} "
                  f"(ML tree: {summ.get('inferred_direction_ml', '?')}) | archaeal clades "
                  f"{summ.get('n_archaeal_clades', '?')} (supported "
                  f"{summ.get('n_archaeal_clades_supported', '?')}) | transitions "
                  f"{summ.get('n_crossdomain_transitions', '?')}")
    txt = open(args.tree).read().strip() if os.path.exists(args.tree) else ""
    if not txt or txt in ("();", ";") or txt.count(",") < 2:
        reason = next((ln.split("\t")[2] for ln in open(args.summary)
                       if ln.startswith("status\t") and ln.count("\t") >= 2), "no tree")
        placeholder(args.out_prefix, f"{title}\n({reason.strip()})", formats)
        utils.log("plot_crossdomain_tree: trivial tree; placeholder.")
        return

    dom, src = load_tips(args.tips)
    parent, children, label, length = parse_newick(txt)
    leaves = [v for v in range(len(parent)) if not children[v]]
    missing = [label[v] for v in leaves if label[v] not in dom]
    if missing:
        raise SystemExit(f"plot_crossdomain_tree: {len(missing)} tree tips absent from "
                         f"--tips (e.g. {missing[:3]})")
    x, y, order = layout(parent, children, length)
    comp = {}
    for v in reversed(order):
        if not children[v]:
            comp[v] = {tip_class(label[v], dom, src)}
        else:
            comp[v] = set().union(*(comp[c] for c in children[v]))

    def colour(v):
        c = comp[v]
        if c == {"Reference"}:
            return COLORS["Reference"]
        if c <= {"Bacteria", "Reference"}:
            return COLORS["Bacteria"]
        if c in ({"Archaea"}, {"Paralog"}):
            return COLORS[next(iter(c))]
        return COLORS["mixed"]

    n_tips = len(leaves)
    height = max(3.0, min(0.12 * n_tips + 1.0, 24.0))
    fig, ax = plt.subplots(figsize=(7.2, height))
    segs, cols = [], []
    for v in order:
        if children[v]:
            ys = [y[c] for c in children[v]]
            segs.append([(x[v], min(ys)), (x[v], max(ys))])
            cols.append(colour(v))
            for c in children[v]:
                segs.append([(x[v], y[c]), (x[c], y[c])])
                cols.append(colour(c))
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=0.6))
    for cls in ("Archaea", "Bacteria", "Reference", "Paralog"):
        vs = [v for v in leaves if tip_class(label[v], dom, src) == cls]
        if vs:
            ax.scatter([x[v] for v in vs], [y[v] for v in vs], s=6 if n_tips > 200 else 12,
                       color=COLORS[cls], zorder=3, linewidths=0)
    sup_x, sup_y = [], []
    for v in order:
        if children[v] and parent[v] is not None:
            s = parse_support(label[v])
            if s and (s[0] is None or s[0] >= args.min_alrt) and (s[1] is None or s[1] >= args.min_ufboot):
                sup_x.append(x[v])
                sup_y.append(y[v])
    if sup_x:
        ax.scatter(sup_x, sup_y, s=4, color="black", zorder=4, linewidths=0)
    if n_tips <= args.max_labels:
        xmax = max(x)
        for v in leaves:
            ax.text(x[v] + 0.01 * xmax, y[v], label[v], fontsize=4, va="center",
                    color=COLORS[tip_class(label[v], dom, src)])
    ax.set_ylim(-1, n_tips)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("substitutions / site")
    ax.set_title(title, fontsize=6)
    present = sorted({tip_class(label[v], dom, src) for v in leaves},
                     key=["Archaea", "Bacteria", "Reference", "Paralog"].index)
    names = {"Archaea": "Archaea", "Bacteria": "Bacteria", "Reference": "Prometheus reference",
             "Paralog": "paralog outgroup"}
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=COLORS[k],
                      markersize=5, label=names[k]) for k in present]
    handles.append(Line2D([0], [0], color=COLORS["mixed"], lw=1, label="mixed lineage"))
    handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="black", markersize=3,
                          label=f"SH-aLRT>={args.min_alrt:g} & UFBoot>={args.min_ufboot:g}"))
    ax.legend(handles=handles, fontsize=5, loc="lower left", frameon=False)
    save_fig(fig, args.out_prefix, formats)
    utils.log(f"plot_crossdomain_tree: drew {n_tips} tips ({len(sup_x)} supported nodes).")


if __name__ == "__main__":
    main()
