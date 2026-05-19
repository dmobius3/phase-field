#!/usr/bin/env python3
"""Phase 1 dry-run harness.

Generates synthetic SPARC-format data (a master .mrt table and
*_rotmod.dat rotation-curve files), runs the real run_pipeline.py and
make_figures.py against it, and checks they complete and produce sane
output. This exercises the loaders, the MRT parser, the 27-cell loop,
and the JSON/figure writes before any real SPARC data is downloaded, so
Phase 1 can be a clean single run.

The synthetic data is NOT SPARC. It is written into data/ (gitignored),
used, and removed on success. The harness refuses to run if data/
already holds a master table, so real data is never clobbered.

    python scripts/dryrun.py            run, then clean up
    python scripts/dryrun.py --keep     run, leave the synthetic data
"""

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
from astropy.table import Table  # noqa: E402

from coherence.synthetic import make_galaxy  # noqa: E402

DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
MASTER = DATA / "SPARC_Lelli2016c.mrt"


def _write_rotmod(name, galaxy):
    """Write one synthetic *_rotmod.dat in the 8-column SPARC format."""
    sb = 80.0 * np.exp(-galaxy["rad"] / max(galaxy["R_d"], 0.2))
    lines = ["# Synthetic dry-run rotation curve -- NOT SPARC data",
             "#  Rad  Vobs  errV  Vgas  Vdisk  Vbul  SBdisk  SBbul"]
    for i in range(len(galaxy["rad"])):
        lines.append("%9.4f %9.3f %8.3f %9.3f %9.3f %8.3f %10.3f %10.3f" % (
            galaxy["rad"][i], galaxy["vobs"][i], galaxy["e_vobs"][i],
            galaxy["vgas"][i], galaxy["vdisk"][i], galaxy["vbul"][i],
            sb[i], 0.0))
    (DATA / (name + "_rotmod.dat")).write_text("\n".join(lines) + "\n")


def _specs():
    """(name, kind, shape, curve v_c, master Vflat, quality, inclination)."""
    shapes = ("sharp", "gradual", "bumpy")
    specs = []
    for i, v in enumerate(np.linspace(140.0, 285.0, 22)):  # flat, measurable
        specs.append(("DryFlatM%02d" % i, "flat", shapes[i % 3],
                      float(v), float(v), 1 + i % 2, 45.0 + (i % 5) * 7.0))
    for j, v in enumerate(np.linspace(85.0, 128.0, 6)):     # flat, below cut
        specs.append(("DryFlatL%02d" % j, "flat", shapes[j % 3],
                      float(v), float(v), 1 + j % 2, 40.0 + j * 6.0))
    for k, v in enumerate(np.linspace(100.0, 210.0, 6)):    # rising
        # master Vflat set above the curve max so R_flat is undefined
        specs.append(("DryRise%02d" % k, "rising", "gradual",
                      float(v), float(v) * 1.3, 1 + k % 2, 50.0 + k * 4.0))
    for m, v in enumerate((160.0, 200.0, 240.0)):           # quality 3
        specs.append(("DryQ3_%02d" % m, "flat", "gradual", v, v, 3, 60.0))
    for m, v in enumerate((170.0, 220.0)):                  # low inclination
        specs.append(("DryLowInc%02d" % m, "flat", "gradual", v, v, 1, 22.0))
    return specs


def build_synthetic():
    """Write the synthetic master table and rotation-curve files."""
    DATA.mkdir(exist_ok=True)
    rows = []
    for name, kind, shape, v_curve, vflat, quality, inc in _specs():
        galaxy = make_galaxy(kind, shape=shape, v_c=v_curve)
        _write_rotmod(name, galaxy)
        m_hi = 0.3 * (v_curve / 100.0) ** 3
        l36 = (v_curve / 100.0) ** 4
        rows.append((name, quality, inc, vflat, galaxy["R_d"],
                     galaxy["R_eff"], m_hi, l36))
    # edge case: a master row with no rotation-curve file
    rows.append(("DryGhost", 1, 60.0, 180.0, 2.5, 4.2, 1.0, 3.0))
    # edge case: a rotation-curve file with no master row
    _write_rotmod("DryOrphan", make_galaxy("flat", v_c=150.0))

    cols = list(zip(*rows))
    table = Table()
    table["Galaxy"] = list(cols[0])
    table["Q"] = np.array(cols[1], dtype=int)
    table["Inc"] = np.array(cols[2], dtype=float)
    table["Vflat"] = np.array(cols[3], dtype=float)
    table["Rdisk"] = np.array(cols[4], dtype=float)
    table["Reff"] = np.array(cols[5], dtype=float)
    table["MHI"] = np.array(cols[6], dtype=float)
    table["L[3.6]"] = np.array(cols[7], dtype=float)
    table.write(MASTER, format="mrt", overwrite=True)
    return len(rows)


def _run(script):
    return subprocess.run([sys.executable, str(ROOT / "scripts" / script)],
                          capture_output=True, text=True)


def cleanup():
    for pattern_dir, pat in ((DATA, "*_rotmod.dat"),):
        for p in pattern_dir.glob(pat):
            p.unlink()
    if MASTER.exists():
        MASTER.unlink()
    for directory in (RESULTS, FIGURES):
        for p in directory.iterdir():
            if p.name != ".gitkeep":
                p.unlink()


def main():
    if MASTER.exists():
        print("ABORT: data/SPARC_Lelli2016c.mrt already exists. Refusing to")
        print("run so real data is not clobbered. Clear data/ if this is")
        print("leftover synthetic data.")
        sys.exit(1)

    keep = "--keep" in sys.argv
    ok = True
    try:
        n = build_synthetic()
        print("built %d master rows + 1 orphan rotation-curve file" % n)

        rp = _run("run_pipeline.py")
        print("\n--- run_pipeline.py (exit %d) ---" % rp.returncode)
        print(rp.stdout.strip() or "(no stdout)")
        if rp.returncode != 0:
            print("STDERR:\n" + rp.stderr.strip())
            ok = False

        if ok:
            grid = json.loads((RESULTS / "sensitivity_grid.json").read_text())
            detail = json.loads((RESULTS / "primary_cell.json").read_text())
            cont = grid["prediction_4_contingency"]
            print("\nsensitivity_grid.json: %d cells, all_stable=%s"
                  % (grid["n_cells"], grid["all_stable"]))
            print("closure_identity: %s" % grid["closure_identity"])
            print("contingency: table=%s fisher_p=%.3g auc=%.3f"
                  % (cont["table"], cont["fisher_exact_p"], cont["auc"]))
            print("primary_cell.json: %d galaxy records" % len(detail))
            if grid["n_cells"] != 27:
                print("WARN: expected 27 grid cells")
                ok = False

        if ok:
            mf = _run("make_figures.py")
            print("\n--- make_figures.py (exit %d) ---" % mf.returncode)
            print(mf.stdout.strip() or "(no stdout)")
            if mf.returncode != 0:
                print("STDERR:\n" + mf.stderr.strip())
                ok = False
            else:
                figs = sorted(p.name for p in FIGURES.glob("*.png"))
                print("figures written: %s" % figs)
                if len(figs) != 3:
                    print("WARN: expected 3 figures")
                    ok = False
    finally:
        if ok and not keep:
            cleanup()
            print("\ncleaned up synthetic data and generated outputs")
        else:
            print("\nsynthetic data left in data/ (and outputs) for "
                  "inspection; clear data/ before any real fetch")

    print("\nDRY RUN: %s" % ("PASS" if ok else "FAIL"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
