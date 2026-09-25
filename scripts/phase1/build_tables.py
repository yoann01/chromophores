"""Phase 1: build the homogeneous LUT and the layered 4D table (resumable).

On a farm, run `python build_tables.py layered --slices i j ...` on several
machines against a shared --work directory, then `python build_tables.py merge`.
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from chromophores import homogeneous_lut, table  # noqa: E402

DATA = os.path.join(os.path.dirname(table.__file__), "data")

ap = argparse.ArgumentParser()
ap.add_argument("what", choices=["homogeneous", "layered", "merge", "all"])
ap.add_argument("--work", default=os.path.join(os.path.dirname(__file__), "..", "..", "build", "layered_slices"))
ap.add_argument("--slices", type=int, nargs="*")
ap.add_argument("--photons-layered", type=int, default=30_000)
ap.add_argument("--photons-homogeneous", type=int, default=50_000)
args = ap.parse_args()

t0 = time.time()
if args.what in ("homogeneous", "all"):
    homogeneous_lut.build(os.path.join(DATA, homogeneous_lut.LUT_FILE), args.photons_homogeneous)
    print(f"homogeneous LUT done ({time.time() - t0:.0f} s)", flush=True)
if args.what in ("layered", "all"):
    table.build(args.work, args.photons_layered, args.slices)
    print(f"layered slices done ({time.time() - t0:.0f} s)", flush=True)
if args.what in ("merge", "all"):
    table.merge_slices(args.work, os.path.join(DATA, table.TABLE_FILE))
    print("merged", flush=True)
