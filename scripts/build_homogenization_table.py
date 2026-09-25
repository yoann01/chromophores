"""Precompute the dimensionless two-layer skin table used by chromophores.homogenize."""

import os
import time

from common import FIG_DIR  # noqa: F401  (sets sys.path)

from chromophores import homogenize

t0 = time.time()
table = homogenize.build_table(n_photons=30_000)
path = os.path.join(os.path.dirname(homogenize.__file__), "data", homogenize.TABLE_FILE)
homogenize.save_table(table, path)
print(f"wrote {path} in {time.time() - t0:.0f} s")
