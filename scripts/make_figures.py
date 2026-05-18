#!/usr/bin/env python3
"""Generate the three figures from results/ into figures/ (gitignored).

Phase 2 code: requires results/primary_cell.json from run_pipeline.py.

  figure 1  eta_flat histogram, with the precision window [0.75, 1.25]
            and the coarse gate [0.5, 2.0] both drawn
  figure 2  r_t vs L_f, with the OLS fit and the null-model correlations
  figure 3  ROC curve for the trigger index (AUC and the 0.7 gate), with
            the 2x2 confusion matrix at the T/T_c = 1 boundary
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from coherence.nulls import NULL_MODELS, null_predictor, ols_fit  # noqa: E402
from coherence.nulls import roc_auc, spearman  # noqa: E402

RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
PRECISION_WINDOW = (0.75, 1.25)
COARSE_GATE = (0.5, 2.0)
AUC_MIN = 0.7


def _column(rows, key, mask=None):
    vals = [r[key] for r in rows]
    if mask is not None:
        vals = [v for v, m in zip(vals, mask) if m]
    return np.array([v for v in vals if v is not None], dtype=float)


def figure_eta_flat(rows):
    """Figure 1: eta_flat histogram with both acceptance gates drawn."""
    measurable = [r for r in rows if r["rflat_measurable"]]
    eta = _column(measurable, "eta_flat")
    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.axvspan(*COARSE_GATE, color="0.88", zorder=0,
               label="coarse gate [0.5, 2.0]")
    ax.axvspan(*PRECISION_WINDOW, color="#f0c419", alpha=0.40, zorder=1,
               label="precision window [0.75, 1.25]")
    ax.hist(eta, bins=24, color="#4c72b0", edgecolor="white", zorder=2)
    ax.axvline(1.0, color="k", ls="--", lw=1.0, zorder=3)
    ax.axvline(float(np.median(eta)), color="#c44e52", lw=2.0, zorder=3,
               label=f"median = {np.median(eta):.2f}")
    ax.set_xlabel(r"$\eta_\mathrm{flat} = R_\mathrm{flat}\,/\,L_f$")
    ax.set_ylabel("galaxies")
    ax.set_title(f"Flat-onset ratio  (n = {len(eta)} measurable)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_eta_flat_histogram.png", dpi=150)
    plt.close(fig)


def figure_rt_vs_lf(rows):
    """Figure 2: r_t vs L_f with the OLS fit and null-model correlations."""
    flat = [r for r in rows if r["flat"] and r["r_t"] is not None]
    L_f = _column(flat, "L_f")
    r_t = _column(flat, "r_t")
    fit = ols_fit(L_f, r_t)

    fig, (ax, axn) = plt.subplots(1, 2, figsize=(11.0, 4.6),
                                  gridspec_kw={"width_ratios": [3, 2]})
    ax.scatter(L_f, r_t, s=22, color="#4c72b0", alpha=0.8)
    line = np.array([L_f.min(), L_f.max()])
    ax.plot(line, fit["slope"] * line + fit["intercept"], color="#c44e52",
            lw=2.0, label=f"OLS slope = {fit['slope']:.2f}")
    ax.plot(line, line, color="k", ls="--", lw=1.0, label="slope 1")
    ax.set_xlabel(r"$L_f = v_c^2 / a_0$  [kpc]")
    ax.set_ylabel(r"$r_t$  [kpc]")
    ax.set_title("Transition radius vs coherence scale")
    ax.legend()

    catalog = {k: _column(flat, k) for k in ("v_c", "M_b", "R_d", "R_eff",
                                             "R_last")}
    names = ["L_f"] + list(NULL_MODELS)
    rhos = [spearman(r_t, L_f)[0]]
    for name in NULL_MODELS:
        rhos.append(spearman(r_t, null_predictor(name, catalog))[0])
    colors = ["#c44e52"] + ["#4c72b0"] * len(NULL_MODELS)
    axn.barh(range(len(names)), np.abs(rhos), color=colors)
    axn.set_yticks(range(len(names)))
    axn.set_yticklabels(names)
    axn.invert_yaxis()
    axn.set_xlabel(r"$|\,$Spearman $\rho$ with $r_t\,|$")
    axn.set_title("L_f vs null models")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_rt_vs_Lf.png", dpi=150)
    plt.close(fig)


def _roc_curve(scores, labels):
    order = np.argsort(-np.asarray(scores, dtype=float))
    lab = np.asarray(labels, dtype=bool)[order]
    n_pos, n_neg = lab.sum(), (~lab).sum()
    tpr = np.concatenate([[0.0], np.cumsum(lab) / n_pos])
    fpr = np.concatenate([[0.0], np.cumsum(~lab) / n_neg])
    return fpr, tpr


def figure_trigger_roc(rows):
    """Figure 3: ROC for the trigger index, plus the 2x2 confusion matrix."""
    scores = np.array([r["trigger_ratio"] for r in rows], dtype=float)
    labels = np.array([r["flat"] for r in rows], dtype=bool)
    auc = roc_auc(scores, labels)
    fpr, tpr = _roc_curve(scores, labels)

    fig, (ax, axc) = plt.subplots(1, 2, figsize=(10.5, 4.6),
                                  gridspec_kw={"width_ratios": [3, 2]})
    ax.plot(fpr, tpr, color="#4c72b0", lw=2.0,
            label=f"AUC = {auc:.2f}")
    ax.plot([0, 1], [0, 1], color="k", ls="--", lw=1.0, label="chance")
    ax.set_xlabel("false positive rate")
    ax.set_ylabel("true positive rate")
    ax.set_title(f"Trigger index ROC  (gate: AUC ≥ {AUC_MIN})")
    ax.text(0.55, 0.10, "PASS" if auc >= AUC_MIN else "FAIL",
            color=("#2a9d4a" if auc >= AUC_MIN else "#c44e52"),
            fontsize=13, fontweight="bold")
    ax.legend(loc="lower right")

    # confusion matrix at the T/T_c = 1 decision boundary
    predicted = scores >= 1.0
    cm = np.array([[np.sum(~predicted & ~labels), np.sum(predicted & ~labels)],
                   [np.sum(~predicted & labels), np.sum(predicted & labels)]])
    axc.imshow(cm, cmap="Blues")
    for i in range(2):
        for j in range(2):
            axc.text(j, i, str(cm[i, j]), ha="center", va="center",
                     fontsize=13)
    axc.set_xticks([0, 1])
    axc.set_xticklabels(["pred. rising", "pred. flat"])
    axc.set_yticks([0, 1])
    axc.set_yticklabels(["rising", "flat"])
    axc.set_title(r"Confusion at $T/T_c = 1$")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_trigger_roc.png", dpi=150)
    plt.close(fig)


def main():
    detail = RESULTS / "primary_cell.json"
    if not detail.exists():
        print("No results/primary_cell.json. make_figures.py is Phase 2 code.")
        print("Run scripts/run_pipeline.py first.")
        sys.exit(1)
    rows = json.loads(detail.read_text())
    FIGURES.mkdir(exist_ok=True)
    figure_eta_flat(rows)
    figure_rt_vs_lf(rows)
    figure_trigger_roc(rows)
    print(f"three figures written to {FIGURES}")


if __name__ == "__main__":
    main()
