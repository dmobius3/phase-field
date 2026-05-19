#!/usr/bin/env python3
"""Phase 1: download the SPARC distribution into data/.

This script makes data contact. It refuses to run without the explicit
--confirm-data-contact flag, so that cloning the repository and running
the scripts in order cannot trigger data contact before the
pre-registration (registration/PREREGISTRATION.md) has been read and the
pipeline tagged. The flag is the visible gate between Phase 0 and Phase 1.

SPARC: Lelli, McGaugh & Schombert (2016), AJ 152, 157. The data is not
redistributed in this repository; it is fetched from the upstream site.
"""

import argparse
import sys
import urllib.request
import zipfile
from pathlib import Path

SPARC_BASE = "https://astroweb.cwru.edu/SPARC/"
MASTER_TABLE = "SPARC_Lelli2016c.mrt"
ROTMOD_ARCHIVE = "Rotmod_LTG.zip"

DATA = Path(__file__).resolve().parent.parent / "data"


def _fetch(name, dest):
    url = SPARC_BASE + name
    print(f"  downloading {url}")
    urllib.request.urlretrieve(url, dest)


def main():
    parser = argparse.ArgumentParser(
        description="Download SPARC into data/ (Phase 1 data contact).")
    parser.add_argument(
        "--confirm-data-contact", action="store_true",
        help="required; acknowledges this script makes Phase 1 data contact")
    args = parser.parse_args()

    if not args.confirm_data_contact:
        print("REFUSING TO RUN.")
        print("fetch_data.py makes Phase 1 data contact with SPARC.")
        print("Before running it: read registration/PREREGISTRATION.md and")
        print("confirm the pipeline is tagged v1.0-preregistration.")
        print("Then re-run with:  python scripts/fetch_data.py "
              "--confirm-data-contact")
        sys.exit(1)

    DATA.mkdir(exist_ok=True)
    print(f"fetching SPARC into {DATA}")
    _fetch(MASTER_TABLE, DATA / MASTER_TABLE)

    archive = DATA / ROTMOD_ARCHIVE
    _fetch(ROTMOD_ARCHIVE, archive)
    with zipfile.ZipFile(archive) as zf:
        zf.extractall(DATA)
    archive.unlink()

    n_curves = len(list(DATA.glob("*_rotmod.dat")))
    print(f"done: master table + {n_curves} rotation curves in {DATA}")


if __name__ == "__main__":
    main()
