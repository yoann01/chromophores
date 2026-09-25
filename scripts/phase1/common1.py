"""Shared helpers for phase-1 validation scripts."""

import os
import sys

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", ".."))

import numpy as np  # noqa: E402
from common import FIG_DIR, SERIES, plt, save  # noqa: E402,F401

from chromophores.optics import SkinParams  # noqa: E402

REPORT_DIR = os.path.join(HERE, "..", "..", "docs", "phase1")
os.makedirs(REPORT_DIR, exist_ok=True)

# Parameter ranges of the phase-1 model (aligned with Aliaga & Jarabo 2026 where possible).
RANGES = {
    "melanin_fraction": (0.001, 0.5, True),
    "eumelanin_ratio": (0.0, 1.0, False),
    "blood_fraction": (0.002, 0.07, True),
    "oxygen_saturation": (0.5, 1.0, False),
    "epidermis_thickness_cm": (0.002, 0.035, True),
    "scattering_scale": (0.7, 1.4, True),
}


def random_skin(rng):
    kw = {}
    for name, (lo, hi, lg) in RANGES.items():
        x = rng.random()
        kw[name] = float(np.exp(np.log(lo) + x * (np.log(hi) - np.log(lo)))) if lg else lo + x * (hi - lo)
    return SkinParams(**kw)


def write_report(name, text):
    path = os.path.join(REPORT_DIR, name)
    with open(path, "w") as f:
        f.write(text)
    print(text)
    print("wrote", os.path.relpath(path))
