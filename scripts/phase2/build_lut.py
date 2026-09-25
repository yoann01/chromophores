"""Phase 2, step 3: build the RGB -> chromophores 3D LUTs (population and face priors)."""

import os
import sys
import time

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from chromophores import inversion as inv, lut3d  # noqa: E402

DATA = os.path.join(HERE, "..", "..", "chromophores", "data")

if __name__ == "__main__":
    priors = inv.load_priors()
    n = int(os.environ.get("LUT_N", "33"))
    suffix = "" if n == 33 else f"_{n}"
    for label, fname in (("population", f"rgb_lut{suffix}_v1.npz"), ("region:face", f"rgb_lut_face{suffix}_v1.npz")):
        t0 = time.time()
        lut3d.build(priors[label], os.path.join(DATA, fname), n=n)
        print(f"{label}: {time.time() - t0:.0f} s", flush=True)
