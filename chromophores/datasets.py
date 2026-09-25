"""Loaders for measured skin reflectance datasets (docs/phase1/NIST_ISSA.md)."""

from __future__ import annotations

import os

import numpy as np

REPO = os.path.join(os.path.dirname(__file__), "..")
ISSA_XLSX = os.path.join(REPO, "data", "ISSA", "ISSA_17_Jan_2025_Yan_Lu.xlsx")

ISSA_LOCATIONS = {1: "Back of hand", 2: "Cheek", 3: "Cheek bone", 4: "Chin", 5: "Ear lobe", 6: "Forehead", 7: "Inner arm",
                  8: "Neck", 9: "Nose tip", 10: "Outer arm", 11: "Palm", 12: "Ring finger"}
ISSA_FACE = {2, 3, 4, 6, 9}


def load_issa(path=ISSA_XLSX):
    """International Skin Spectra Archive (ISSA, Leeds Skin Database).

    Returns (meta: pandas.DataFrame, wavelengths [nm], reflectance [0-1], shape (n, n_lambda)).
    Reflectance is converted from percent; NaN where a record does not cover a wavelength.
    """
    import pandas as pd

    raw = pd.read_excel(path, sheet_name="ISSA", header=None)
    head = raw.iloc[11]
    lam_cols = [c for c in raw.columns if isinstance(head[c], (int, float, np.number)) and not pd.isna(head[c]) and c >= 13]
    lam = head[lam_cols].to_numpy(float)
    body = raw.iloc[12:]
    body = body[pd.to_numeric(body[0], errors="coerce").notna()]
    meta = body.iloc[:, :12].copy()
    meta.columns = ["record", "origin", "subject", "ethnicity", "gender", "age_group", "location", "instrument", "spin",
                    "start_nm", "end_nm", "step_nm"]
    meta = meta.reset_index(drop=True)
    R = body[lam_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float) / 100.0
    return meta, lam, R


NIST_TXT = os.path.join(REPO, "data", "nist", "1832_Data_JResNIST_skinrefl v3.txt")


def load_nist(path=NIST_TXT, average=True):
    """NIST Reference Data Set of Human Skin Reflectance (Cooksey et al. 2017, doi:10.18434/M38597).

    100 subjects, inner right forearm, directional-hemispherical reflectance factor (integrating
    sphere, specular included), 250-2500 nm / 3 nm. Returns (wavelengths [nm], reflectance
    (n_subjects, n_lambda)) using the per-subject average of the 3 repeats (or all repeats).
    """
    import csv

    with open(path, encoding="latin-1") as f:
        rows = list(csv.reader(f))
    labels = rows[7]
    data = np.array([[float(v) if v else np.nan for v in r[: len(labels)]] for r in rows[8:] if r and r[0]])
    lam = data[:, 0]
    cols = [i for i, lab in enumerate(labels) if (lab == "Average") == average and lab in ("R1", "R2", "R3", "Average")]
    return lam, data[:, cols].T
