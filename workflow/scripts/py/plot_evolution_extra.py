"""Extra evolution figures (added alongside plot_selection.py / the R evolution plots).

Outputs (written to --outdir):
  dnds_along_protein.{png,svg}      per-site dN/dS (FEL/MEME) along the alignment with a
                                    domain track and selected sites flagged -- shows WHERE
                                    selection acts (diversifying signal on the nuclease/
                                    REase modules = the host-virus conflict signature).
  reconciliation_events.{png,svg}   D/T/L/S event counts per reconciliation tool + the
                                    vertical (speciation/(speciation+transfer)) fraction.
  gain_enrichment_null.{png,svg}    methanogen-gain-enrichment: the phylogenetically
                                    structured null distribution with the observed value
                                    and one-sided p (honest view of the enrichment test).

These foreground the conflict/HGT/enrichment evidence with their uncertainty rather than
as single numbers.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

apply_style()
# Domain-track colours chosen to NOT collide with the selection-class colours
# (positive=#D55E00 vermillion, negative=#0072B2 blue, neutral=#999999 grey).
DOMAIN_COLORS = {"DUF4011": "#E69F00", "DUF3320": "#56B4E9", "DNA2": "#009E73",
                 "nuclease": "#CC79A7", "REase": "#882255", "NA": "#DDDDDD"}


def _load(path):
    if not path or not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return None
    return None if df.empty else df


def _load_null(path):
    """Load the enrichment-null file WITHOUT pandas' default NaN coercion, so a 'kind'
    label is never silently turned into NaN (pandas reads 'null'/'NA'/'None' as NaN)."""
    if not path or not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[])
    except Exception:
        return None
    return None if df.empty else df


# --------------------------------------------------------------------------- #
def plot_dnds_along_protein(sites, outdir):
    prefix = os.path.join(outdir, "dnds_along_protein")
    if sites is None or "alignment_col" not in sites.columns:
        placeholder(prefix, "no data: site_table.tsv missing/empty")
        return
    df = sites.copy()
    df["x"] = pd.to_numeric(df["alignment_col"], errors="coerce")
    df["stat"] = pd.to_numeric(df.get("beta_alpha_or_dnds", np.nan), errors="coerce")
    df = df.dropna(subset=["x"])
    if df.empty:
        placeholder(prefix, "no data: no usable site rows")
        return
    cls = df.get("classification", pd.Series("neutral", index=df.index)).astype(str)
    cmap = {"positive": OKABE_ITO["vermillion"], "negative": OKABE_ITO["blue"], "neutral": OKABE_ITO["grey"]}
    colors = cls.map(lambda c: cmap.get(c, OKABE_ITO["grey"]))
    sizes = np.where(cls.values == "positive", 16, 6)

    fig, (ax, axd) = plt.subplots(2, 1, figsize=(7.0, 3.0), sharex=True,
                                  gridspec_kw={"height_ratios": [6, 1], "hspace": 0.06})
    yvals = df["stat"].fillna(1.0).clip(lower=1e-3)
    ax.scatter(df["x"].values, yvals.values, c=colors.values, s=sizes, edgecolor="none", zorder=3)
    ax.axhline(1.0, color="black", linewidth=0.6, linestyle="--")
    ax.set_yscale("log")
    ax.set_ylabel("dN/dS (β/α)")
    ax.set_title("Per-site selection along the Prometheus alignment")
    # domain track
    if "domain" in df.columns:
        dom = df.sort_values("x")[["x", "domain"]]
        for _, r in dom.iterrows():
            axd.add_patch(plt.Rectangle((r["x"] - 0.5, 0), 1, 1,
                          color=DOMAIN_COLORS.get(str(r["domain"]), DOMAIN_COLORS["NA"]), linewidth=0))
        axd.set_xlim(df["x"].min() - 1, df["x"].max() + 1)
    axd.set_yticks([]); axd.set_ylim(0, 1); axd.set_xlabel("alignment column")
    handles = [mpatches.Patch(color=cmap["positive"], label="positive"),
               mpatches.Patch(color=cmap["negative"], label="negative"),
               mpatches.Patch(color=cmap["neutral"], label="neutral")]
    dom_handles = [mpatches.Patch(color=v, label=k) for k, v in DOMAIN_COLORS.items() if k != "NA"]
    ax.legend(handles=handles + dom_handles, frameon=False, fontsize=5, ncol=2, loc="upper right")
    save_fig(fig, prefix)


# --------------------------------------------------------------------------- #
def plot_reconciliation_events(rec, outdir):
    prefix = os.path.join(outdir, "reconciliation_events")
    if rec is None or not {"tool", "event", "count"}.issubset(rec.columns):
        placeholder(prefix, "no data: reconciliation_events.tsv missing/empty")
        return
    df = rec.copy()
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0)
    events = ["Duplication", "Transfer", "Loss", "Speciation", "Origination"]
    df = df[df["event"].isin(events)]
    tools = sorted(df["tool"].dropna().unique())
    if df.empty or not tools:
        placeholder(prefix, "no data: no recognised reconciliation events")
        return
    x = np.arange(len(events)); w = 0.8 / max(1, len(tools))
    fig, ax = plt.subplots(figsize=(5.0, 2.8))
    palette = [OKABE_ITO["blue"], OKABE_ITO["orange"], OKABE_ITO["green"]]
    for i, tool in enumerate(tools):
        sub = df[df["tool"] == tool].set_index("event")["count"]
        vals = [sub.get(e, 0) for e in events]
        ax.bar(x + i * w, vals, width=w, label=tool, color=palette[i % len(palette)], alpha=0.85)
        # vertical fraction annotation
        frac = df[(df["tool"] == tool) & (df["event"] == "Transfer")].get("fraction")
        fr = ""
        try:
            fv = df[(df["tool"] == tool)]["fraction"].dropna().astype(float)
            if len(fv):
                fr = f"vertical frac={fv.iloc[0]:.2f}"
        except Exception:
            fr = ""
        if fr:
            ax.text(0.02, 0.95 - 0.07 * i, f"{tool}: {fr}", transform=ax.transAxes, fontsize=5)
    ax.set_xticks(x + (len(tools) - 1) * w / 2); ax.set_xticklabels(events, rotation=30, ha="right", fontsize=6)
    ax.set_ylabel("event count"); ax.set_title("DTL reconciliation events")
    ax.legend(frameon=False, fontsize=6)
    save_fig(fig, prefix)


# --------------------------------------------------------------------------- #
def plot_gain_enrichment(gl, null, outdir):
    prefix = os.path.join(outdir, "gain_enrichment_null")
    obs = expp = p = np.nan
    if gl is not None and {"metric", "value"}.issubset(gl.columns):
        vals = dict(zip(gl["metric"].astype(str), pd.to_numeric(gl["value"], errors="coerce")))
        obs = vals.get("gains_on_methanogen_branches", np.nan)
        expp = vals.get("expected_gains_methanogen", np.nan)
        p = vals.get("enrichment_p", np.nan)
    null_vals = None
    if null is not None and "value" in null.columns:
        kind = null.get("kind", pd.Series("sim", index=null.index)).astype(str)
        is_null = kind.isin(["sim", "null", "permutation"])
        null_vals = pd.to_numeric(null.loc[is_null, "value"], errors="coerce").dropna().values
        if "observed" in set(kind) and not np.isfinite(obs):
            obs = pd.to_numeric(null.loc[kind == "observed", "value"], errors="coerce").dropna().mean()
    if (null_vals is None or len(null_vals) == 0) and not np.isfinite(obs):
        placeholder(prefix, "no data: ancestral enrichment not available")
        return

    fig, ax = plt.subplots(figsize=(4.4, 2.8))
    if null_vals is not None and len(null_vals):
        ax.hist(null_vals, bins=min(30, max(8, len(null_vals) // 10)), color=OKABE_ITO["grey"],
                alpha=0.7, label="null (Mk-simulated)")
    elif np.isfinite(expp):
        ax.axvline(expp, color=OKABE_ITO["grey"], linewidth=1.2, label=f"null mean={expp:.2f}")
    if np.isfinite(obs):
        ax.axvline(obs, color=OKABE_ITO["vermillion"], linewidth=1.6, label=f"observed={obs:.2f}")
    ax.set_xlabel("gains on methanogen branches")
    ax.set_ylabel("null simulations" if (null_vals is not None and len(null_vals)) else "")
    title = "Methanogen enrichment of Prometheus gains"
    if np.isfinite(p):
        title += f"  (one-sided p={p:.3g})"
    ax.set_title(title, fontsize=7)
    ax.legend(frameon=False, fontsize=6)
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
    plot_dnds_along_protein(_load(args.sites), args.outdir)
    plot_reconciliation_events(_load(args.reconciliation), args.outdir)
    plot_gain_enrichment(_load(args.gains), _load_null(args.gains_null), args.outdir)


if __name__ == "__main__":
    main()
