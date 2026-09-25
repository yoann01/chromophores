"""Phase 1: fit K=1 and K=2 mixtures on every cell of the layered table."""

import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from chromophores import homogeneous_lut, mixture, table  # noqa: E402

DATA = os.path.join(os.path.dirname(table.__file__), "data")

if __name__ == "__main__":
    t0 = time.time()
    mixture.fit_table(
        os.path.join(DATA, table.TABLE_FILE),
        os.path.join(DATA, homogeneous_lut.LUT_FILE),
        os.path.join(DATA, mixture.FIT_FILE),
        processes=os.cpu_count(),
    )
    print(f"mixture fit done in {time.time() - t0:.0f} s")
