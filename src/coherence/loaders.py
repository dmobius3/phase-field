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

    The master table is a CDS machine-readable table; astropy's MRT
    reader handles the header. Returns a pandas DataFrame.
    """
    from astropy.io import ascii as astropy_ascii

    table = astropy_ascii.read(str(path), format="mrt")
    return table.to_pandas()


def rotation_curve_files(data_dir):
    """All ``*_rotmod.dat`` paths under ``data_dir``, sorted by name."""
    return sorted(Path(data_dir).glob("*_rotmod.dat"))
