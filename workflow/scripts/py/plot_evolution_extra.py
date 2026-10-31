"""Extra evolution figures (rule plot_evolution_extra in plots_extra.smk; CLI frozen).

Outputs (written to --outdir):
  dnds_along_protein.{png,svg}      per-site FEL beta/alpha along the codon alignment with
                                    the domain track (domains taken from the site table, in
                                    N->C order) and selected sites flagged. Sites with alpha=0
                                    and beta>0 (unbounded ratio) are drawn as triangles on the
                                    top edge; sites with no information (0/0) are not drawn.
  reconciliation_events.{png,svg}   D/T/L/S event counts per reconciliation tool + the
                                    vertical fraction (S+SL)/(S+SL+T+TL) read from the
                                    Speciation row; tools that did not estimate are labelled.
  gain_enrichment_null.{png,svg}    methanogen-gain enrichment for the PRIMARY tier set (first
                                    tier_set in the file): the structured null, the observed
                                    value and the one-sided p.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO  # noqa: E402
import argparse  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402

apply_style()
CLASS_COLORS = {"positive": OKABE_ITO["vermillion"], "negative": OKABE_ITO["blue"],
                "neutral": OKABE_ITO["grey"]}
# domain-track palette avoids the class colours above
DOMAIN_PALETTE = [OKABE_ITO["orange"], OKABE_ITO["green"], "#882255", OKABE_ITO["skyblue"],
                  OKABE_ITO["purple"], OKABE_ITO["yellow"], "#44AA99", "#332288"]
EVENTS = ["Duplication", "Transfer", "Loss", "Speciation", "Origination"]


def _load(path, **kw):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    df = pd.read_csv(path, sep="\t", **kw)
    return None if df.empty else df


def _load_null(path):
    """Keep literal labels (pandas would turn 'null'/'NA' into NaN)."""
    return _load(path, keep_default_na=False, na_values=[])


def _primary(df):
    if df is None or "tier_set" not in df.columns:
        return df
    first = df["tier_set"].astype(str).iloc[0]
    return df[df["tier_set"].astype(str) == first]


# --------------------------------------------------------------------------- #
def plot_dnds_along_protein(sites, outdir):
    prefix = os.path.join(outdir, "dnds_along_protein")
    if sites is None or "alignment_col" not in sites.columns:
        placeholder(prefix, "no data: site_table.tsv missing/empty (selection not run)")
        return
    df = sites.copy()
    df["x"] = pd.to_numeric(df["alignment_col"], errors="coerce")
    raw = df.get("beta_alpha_or_dnds", pd.Series(np.nan, index=df.index)).astype(str).str.strip()
    df["stat"] = pd.to_numeric(raw.replace({"inf": "inf"}), errors="coerce")
    df = df.dropna(subset=["x"]).sort_values("x")
    if df.empty:
        placeholder(prefix, "no data: no usable site rows")
        return
    cls = df.get("classification", pd.Series("neutral", index=df.index)).fillna("neutral").astype(str)
    finite = np.isfinite(df["stat"]) & (df["stat"] > 0)
    inf = np.isinf(df["stat"])
    fig, (ax, axd) = plt.subplots(2, 1, figsize=(7.0, 3.0), sharex=True,
                                  gridspec_kw={"height_ratios": [6, 1], "hspace": 0.06})
    ax.scatter(df.loc[finite, "x"], df.loc[finite, "stat"], c=cls[finite].map(CLASS_COLORS).fillna(OKABE_ITO["grey"]),
               s=np.where(cls[finite] == "positive", 14, 5), edgecolor="none", zorder=3)
    ax.set_yscale("log")
    if finite.any():
        top = max(10.0, float(df.loc[finite, "stat"].max()) * 1.5)
    else:
        top = 10.0
    ax.set_ylim(min(0.01, float(df.loc[finite, "stat"].min()) / 1.5) if finite.any() else 0.01, top)
    if inf.any():
        ax.scatter(df.loc[inf, "x"], np.full(int(inf.sum()), top / 1.15), marker="^",
                   c=cls[inf].map(CLASS_COLORS).fillna(OKABE_ITO["grey"]), s=10, edgecolor="none",
                   zorder=3, clip_on=False)
    ax.axhline(1.0, color="black", linewidth=0.6, linestyle="--")
    ax.set_ylabel("FEL β/α (log)")
    ax.set_title(f"Per-site selection along the Prometheus codon alignment "
                 f"({int((cls == 'positive').sum())} positive, {int((cls == 'negative').sum())} negative; "
                 f"triangles: α=0, β>0)", fontsize=7)
    dom_handles = []
    if "domain" in df.columns:
        doms = [d for d in dict.fromkeys(df["domain"].fillna("NA").astype(str)) if d != "NA"]
        cmap = {d: DOMAIN_PALETTE[i % len(DOMAIN_PALETTE)] for i, d in enumerate(doms)}
        cmap["NA"] = "#EEEEEE"
        for x, d in zip(df["x"], df["domain"].fillna("NA").astype(str)):
            axd.add_patch(plt.Rectangle((x - 0.5, 0), 1, 1, color=cmap.get(d, "#EEEEEE"), linewidth=0))
        dom_handles = [mpatches.Patch(color=cmap[d], label=d) for d in doms]
    axd.set_xlim(df["x"].min() - 1, df["x"].max() + 1)
    axd.set_yticks([]); axd.set_ylim(0, 1); axd.set_xlabel("codon alignment column")
    handles = [mpatches.Patch(color=CLASS_COLORS[k], label=k) for k in ("positive", "negative", "neutral")]
    ax.legend(handles=handles + dom_handles, frameon=False, fontsize=5, ncol=2, loc="upper right")
    save_fig(fig, prefix)


# --------------------------------------------------------------------------- #
def plot_reconciliation_events(rec, outdir):
    prefix = os.path.join(outdir, "reconciliation_events")
    if rec is None or not {"tool", "event", "count"}.issubset(rec.columns):
        placeholder(prefix, "no data: reconciliation_events.tsv missing/empty")
        return
    df = rec[rec["event"].isin(EVENTS)].copy()
    df["count"] = pd.to_numeric(df["count"], errors="coerce")
    notes = df.get("note", pd.Series("", index=df.index)).fillna("").astype(str)
    est = ~notes.str.contains(r"disabled|too few|not run|trivial", case=False, regex=True)
    tools = list(dict.fromkeys(df["tool"].astype(str)))
    if df.empty or not tools:
        placeholder(prefix, "no data: no recognised reconciliation events")
        return
    x = np.arange(len(EVENTS)); w = 0.8 / max(1, len(tools))
    fig, ax = plt.subplots(figsize=(5.0, 2.8))
    palette = [OKABE_ITO["blue"], OKABE_ITO["orange"], OKABE_ITO["green"]]
    for i, tool in enumerate(tools):
        sub = df[(df["tool"].astype(str) == tool)]
        ok = bool(est[sub.index].all())
        cnt = sub.set_index("event")["count"]
        vals = [float(cnt.get(e, 0) or 0) if ok else 0 for e in EVENTS]
        ax.bar(x + i * w, vals, width=w, label=tool if ok else f"{tool} (not estimated)",
               color=palette[i % len(palette)], alpha=0.85 if ok else 0.25)
        txt = "not estimated: " + notes[sub.index[0]][:60] if not ok else ""
        if ok:
            sp = sub[sub["event"] == "Speciation"]
            fv = pd.to_numeric(sp.get("fraction", pd.Series(dtype=float)), errors="coerce").dropna()
            if len(fv):
                txt = f"vertical fraction (S+SL)/(S+SL+T+TL) = {fv.iloc[0]:.2f}"
        if txt:
            ax.text(0.02, 0.95 - 0.08 * i, f"{tool}: {txt}", transform=ax.transAxes, fontsize=5, va="top")
    ax.set_xticks(x + (len(tools) - 1) * w / 2)
    ax.set_xticklabels(EVENTS, rotation=30, ha="right", fontsize=6)
    ax.set_ylabel("event count"); ax.set_title("DTL reconciliation events")
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    save_fig(fig, prefix)


# --------------------------------------------------------------------------- #
def plot_gain_enrichment(gl, null, outdir):
    prefix = os.path.join(outdir, "gain_enrichment_null")
    gl, null = _primary(gl), _primary(null)
    obs = expp = p = np.nan
    note = ""
    tier = ""
    if gl is not None and {"metric", "value"}.issubset(gl.columns):
        vals = dict(zip(gl["metric"].astype(str), pd.to_numeric(gl["value"], errors="coerce")))
        obs = vals.get("gains_on_methanogen_branches", np.nan)
        expp = vals.get("expected_gains_methanogen", np.nan)
        p = vals.get("enrichment_p", np.nan)
        if "note" in gl.columns:
            n = gl.loc[gl["metric"].astype(str) == "enrichment_p", "note"]
            note = str(n.iloc[0]) if len(n) else ""
        if "tier_set" in gl.columns:
            tier = str(gl["tier_set"].iloc[0])
    null_vals = None
    if null is not None and {"kind", "value"} <= set(null.columns):
        kind = null["kind"].astype(str)
        null_vals = pd.to_numeric(null.loc[kind.isin(["sim", "null", "permutation"]), "value"],
                                  errors="coerce").dropna().values
        if not np.isfinite(obs) and (kind == "observed").any():
            obs = pd.to_numeric(null.loc[kind == "observed", "value"], errors="coerce").dropna().mean()
    if (null_vals is None or len(null_vals) == 0) and not np.isfinite(obs):
        placeholder(prefix, "no data: methanogen-gain enrichment not computed" + (f" ({note[:80]})" if note else ""))
        return
    fig, ax = plt.subplots(figsize=(4.4, 2.8))
    if null_vals is not None and len(null_vals):
        cond = "prevalence-matched" if "prevalence-matched" in note else "unconditional"
        ax.hist(null_vals, bins=min(30, max(8, len(null_vals) // 10)), color=OKABE_ITO["grey"],
                alpha=0.7, label=f"structured null ({cond} Mk), n={len(null_vals)}")
    if np.isfinite(expp):
        ax.axvline(expp, color=OKABE_ITO["black"], linewidth=0.8, linestyle=":", label=f"null mean={expp:.2f}")
    if np.isfinite(obs):
        ax.axvline(obs, color=OKABE_ITO["vermillion"], linewidth=1.6, label=f"observed={obs:.2f}")
    ax.set_xlabel("expected gains on all-methanogen branches")
    ax.set_ylabel("null simulations")
    title = "Methanogen enrichment of Prometheus gains"
    if tier:
        title += f" (tiers {tier})"
    if np.isfinite(p):
        title += f"\none-sided p={p:.3g}"
    ax.set_title(title, fontsize=7)
    ax.legend(frameon=False, fontsize=5.5)
    save_fig(fig, prefix)


def main():
    ap = argparse.ArgumentParser(description="Extra evolution figures")
    ap.add_argument("--sites", required=True)
    ap.add_argument("--reconciliation", required=True)
    ap.add_argument("--gains", required=True)
    ap.add_argument("--gains-null", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    plot_dnds_along_protein(_load(args.sites, keep_default_na=False, na_values=[""]), args.outdir)
    plot_reconciliation_events(_load(args.reconciliation, keep_default_na=False, na_values=[""]), args.outdir)
    plot_gain_enrichment(_load(args.gains, keep_default_na=False, na_values=["", "NA"]),
                         _load_null(args.gains_null), args.outdir)


if __name__ == "__main__":
    main()
