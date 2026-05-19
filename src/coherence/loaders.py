"""Parsers for the SPARC distribution (Lelli, McGaugh & Schombert 2016).

SPARC ships a machine-readable master table and one rotation-curve file
per galaxy. These parsers are frozen with the rest of the pipeline and
run for the first time in Phase 1; nothing here is exercised in Phase 0.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROTMOD_COLUMNS = ["rad", "vobs", "e_vobs", "vgas", "vdisk", "vbul",
                  "sb_disk", "sb_bul"]

# SPARC Table1 (SPARC_Lelli2016c.mrt) columns, in file order.
MASTER_COLUMNS = ["Galaxy", "T", "D", "e_D", "f_D", "Inc", "e_Inc",
                  "L[3.6]", "e_L[3.6]", "Reff", "SBeff", "Rdisk", "SBdisk",
                  "MHI", "RHI", "Vflat", "e_Vflat", "Q", "Ref"]


def load_rotation_curve(path):
    """Parse a SPARC ``*_rotmod.dat`` file.

    Columns: rad (kpc); vobs, e_vobs, vgas, vdisk, vbul (km/s); sb_disk,
    sb_bul (Lsun/pc^2). Comment lines beginning with ``#`` are skipped.
    Returns a DataFrame sorted by radius.
    """
    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None,
                     names=ROTMOD_COLUMNS)
    return df.sort_values("rad").reset_index(drop=True)


def load_master(path):
    """Parse the SPARC master table ``SPARC_Lelli2016c.mrt``.

    The file is an MRT with a byte-by-byte header, but its data rows are
    whitespace-delimited and do not align byte-for-byte with the declared
    column widths, so astropy's strict fixed-width MRT reader fails on it.
    The data section, after the last dashed separator line, is read as
    whitespace-delimited columns in the documented Table1 order
    (``MASTER_COLUMNS``); the trailing reference field absorbs any extra
    tokens. Returns a pandas DataFrame of strings; callers convert the
    columns they use.
    """
    lines = Path(path).read_text().splitlines()
    seps = [i for i, ln in enumerate(lines)
            if ln.strip() and set(ln.strip()) == {"-"}]
    start = seps[-1] + 1 if seps else 0
    n_fixed = len(MASTER_COLUMNS) - 1  # all columns but the trailing Ref
    rows = []
    for line in lines[start:]:
        parts = line.split()
        if len(parts) < len(MASTER_COLUMNS):
            continue
        rows.append(parts[:n_fixed] + [" ".join(parts[n_fixed:])])
    return pd.DataFrame(rows, columns=MASTER_COLUMNS)


def rotation_curve_files(data_dir):
    """All ``*_rotmod.dat`` paths under ``data_dir``, sorted by name."""
    return sorted(Path(data_dir).glob("*_rotmod.dat"))
