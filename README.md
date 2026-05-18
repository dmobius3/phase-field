# phase-field

A pre-registered SPARC test of the Mode Identity Theory coherence scale.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20271702.svg)](https://doi.org/10.5281/zenodo.20271702)

## What this is

Mode Identity Theory (MIT) proposes a galactic coherence scale

    L_f = v_c^2 / a_0

inside which a binary phase field operates. This repository tests one
empirical question against the SPARC sample of 175 disk galaxies:

> Does `v_c^2/a_0` behave like a real coherence radius, after controlling
> for ordinary galactic size scaling?

The theoretical motivation lives in the companion working note
([mode-identity-theory](https://github.com/dmobius3/mode-identity-theory),
`files/working/files/sparc-phase-field.md`). This repository is the
analysis: the pipeline, the frozen pre-registration, and the figures.

## Pre-registration protocol

This is a locked-pipeline blind analysis of pre-existing data. SPARC has
been public since 2016, so the analysis runs under a fixed discipline:

| Phase | Step |
|---|---|
| 0 | Build the full pipeline and synthetic rotation-curve tests. Compute the acceptance scatter bound `sigma_pred` from a representative SPARC error budget. Freeze `registration/PREREGISTRATION.md`. No SPARC data is downloaded. |
| 0 | Tag `v1.0-preregistration` and archive that tag on Zenodo. |
| 1 | Run `scripts/fetch_data.py`, then run the pipeline once. |
| 2 | Generate the figures, the failure-mode table, and the writeup. |

The frozen tag fixes the algorithm thresholds before data contact. It
does not claim the data was unseen. The primary defense is that
`L_f = v_c^2/a_0` and the slope-near-one predictions are parameter-free
consequences of the MIT topology: nothing in the pipeline is fitted. The
genuine pre-registered test against unseen data is Euclid DR1. See
`registration/PREREGISTRATION.md`.

## Layout

    src/coherence/   importable pipeline modules
    scripts/         CLI entry points (fetch_data, run_pipeline, make_figures)
    tests/           synthetic rotation-curve tests
    registration/    frozen pre-registration and acceptance bounds
    data/            SPARC data, fetched by fetch_data.py, not committed
    results/         generated
    figures/         generated

## Reproduce

    python -m venv venv && source venv/bin/activate
    pip install -r requirements.txt
    pytest                          # synthetic-curve tests (no data needed)
    python scripts/fetch_data.py --confirm-data-contact   # Phase 1: download SPARC
    python scripts/run_pipeline.py
    python scripts/make_figures.py

The `--confirm-data-contact` flag on `fetch_data.py` is the explicit
Phase 0 / Phase 1 boundary: the script refuses to run without it, so
cloning the repository and running the scripts cannot trigger data
contact by accident.

## Data

SPARC (Spitzer Photometry and Accurate Rotation Curves): Lelli, McGaugh &
Schombert (2016), AJ 152, 157. Source: astroweb.cwru.edu/SPARC. SPARC is
not redistributed in this repository; `scripts/fetch_data.py` downloads it
from the upstream site.

## License and citation

MIT, see `LICENSE`. Citation metadata in `CITATION.cff`.
