"""Phase 2 production pipeline: albedo texture -> chromophore maps (+ uncertainty, mask, residual).

Steps (ADR-0003, docs/phase2):
1. Auto-calibration (exposure + R/B balance) against the measured skin gamut.
2. 3D LUT lookup: per-texel MAP in unbounded parameter space + full Laplace covariance.
3. Non-skin mask: out of the skin gamut, or chromophores implausible under the prior
   (Mahalanobis distance), e.g. hair (extreme melanin) or mouth interior (extreme blood).
4. Uncertainty-guided spatial regularisation of the poorly constrained parameters only. Melanin and
   blood are kept at their per-texel MAP (full detail; their melanin/blood cross-talk is a
   low-frequency model bias that no local smoothing can fix, see docs/phase2/CLOTURE.md). The
   other parameters o are pulled towards a spatially smoothed estimate, conditionally on the kept
   ones:  x_o = (Q_oo + lam I)^-1 (Q_oo x_o + lam blur(x_o)),  Q = C^-1.
5. Forward model on every texel (maps -> spectrum -> RGB) and colour residual: 3 smooth spectral
   coefficients per texel such that maps + residual reproduce the input exactly.
6. EXR export.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import gaussian_filter

from . import inversion as inv
from . import lut3d, table
from .forward import SPECULAR
from .optics import SkinParams

LAM = inv.LAMBDA
_TABLE = {}


def spectra_from_params(P, chunk=40_000):
    """Vectorised forward model: physical parameters (N, 6) -> diffuse albedo spectra (N, n_lambda)."""
    if "t" not in _TABLE:
        _TABLE["t"] = table.LayeredTable(method="linear")
    T = _TABLE["t"]
    out = np.empty((len(P), len(LAM)))
    for i in range(0, len(P), chunk):
        c = P[i : i + chunk]
        p = SkinParams(bilirubin_umol=0.0, **{n: c[:, k : k + 1] for k, n in enumerate(inv.NAMES)})
        tau, ad, de, g, _ = table.dimensionless(LAM, p)
        out[i : i + chunk] = (1 - SPECULAR) * np.exp(T._A(T.coords(tau, ad, de, np.broadcast_to(g, tau.shape))))
    return out


def colour_residual(S, y, M=lut3d.M_SRGB):
    """Batched `inversion.colour_exact`: coefficients c (N, 3) such that M @ (S * (1 + B^T c)) = y."""
    B = np.abs(M) / np.abs(M).max(1, keepdims=True)
    A = np.einsum("ij,nj,kj->nik", M, S, B)  # (N, 3, 3)
    rhs = y - S @ M.T
    return np.linalg.solve(A, rhs[..., None])[..., 0]


def masked_blur(x, mask, sigma):
    """Gaussian blur of each channel restricted to the mask (normalised convolution)."""
    w = gaussian_filter(mask.astype(float), sigma)
    out = np.empty_like(x)
    for k in range(x.shape[-1]):
        out[..., k] = gaussian_filter(np.where(mask, x[..., k], 0.0), sigma) / np.maximum(w, 1e-6)
    return out


@dataclass
class Result:
    gains: np.ndarray
    albedo: np.ndarray  # calibrated input, linear sRGB
    x: np.ndarray  # regularised MAP, unbounded space (H, W, 6)
    x_raw: np.ndarray  # before regularisation
    std: np.ndarray  # posterior std, unbounded space (H, W, 6)
    skin: np.ndarray  # (H, W) bool
    in_gamut: np.ndarray
    mahalanobis: np.ndarray
    rgb_model: np.ndarray  # maps only (no residual)
    residual: np.ndarray  # (H, W, 3) coefficients
    info: dict = field(default_factory=dict)

    @property
    def maps(self):
        return inv.to_physical(self.x)

    def relative_sigma(self):
        """~1-sigma relative uncertainty of each physical map: log(p84/p16)/2."""
        lo, hi = inv.to_physical(self.x - self.std), inv.to_physical(self.x + self.std)
        return np.log(np.maximum(hi, 1e-12) / np.maximum(lo, 1e-12)) / 2


def process_albedo(albedo_lin, lut_path=None, calibrate=True, prior_label="region:face", mask_chi2=22.46, reg_sigma=6.0,
                   reg_lambda=1.0, reg_keep=("melanin_fraction", "blood_fraction"), chunk=250_000):
    """Albedo (H, W, 3) linear sRGB -> Result. mask_chi2 = chi2(6 dof) 99.9 % quantile."""
    H, W, _ = albedo_lin.shape
    gains = lut3d.autocalibrate(albedo_lin)[0] if calibrate else np.ones(3)
    rgb = (albedo_lin * gains).reshape(-1, 3)
    lut = lut3d.RGBToChromophores(lut_path or os.path.join(os.path.dirname(__file__), "data", "rgb_lut_face_65_v1.npz"))
    if lut._cov is None:
        raise RuntimeError("LUT without covariance: run lut3d.add_covariance(path) first")
    x = np.empty((len(rgb), inv.N))
    C = np.empty((len(rgb), inv.N, inv.N))
    ing = np.empty(len(rgb), bool)
    for i in range(0, len(rgb), chunk):
        o = lut(rgb[i : i + chunk])
        x[i : i + chunk], C[i : i + chunk], ing[i : i + chunk] = o["x"], o["cov"], o["in_gamut"]

    prior = inv.load_priors()[prior_label]
    Pinv = np.linalg.inv(prior.cov)
    d = x - prior.mean
    maha = np.einsum("ni,ij,nj->n", d, Pinv, d)
    skin = ing & (maha <= mask_chi2)

    x_img = x.reshape(H, W, inv.N)
    skin_img = skin.reshape(H, W)
    x_reg = x_img.copy()
    if reg_sigma > 0:
        o = np.array([k for k, n in enumerate(inv.NAMES) if n not in reg_keep])
        x_blur = masked_blur(x_img[..., o], skin_img, reg_sigma).reshape(-1, len(o))
        Q = np.linalg.inv(C + 1e-6 * np.eye(inv.N))[:, o[:, None], o[None, :]]
        rhs = np.einsum("nij,nj->ni", Q, x[:, o]) + reg_lambda * x_blur
        xo = np.linalg.solve(Q + reg_lambda * np.eye(len(o)), rhs[..., None])[..., 0]
        xr = x.copy()
        xr[:, o] = xo
        x_reg = np.where(skin[:, None], xr, x).reshape(H, W, inv.N)

    S = spectra_from_params(inv.to_physical(x_reg.reshape(-1, inv.N)))
    rgb_model = S @ lut3d.M_SRGB.T
    res = colour_residual(S, rgb)
    std = np.sqrt(np.maximum(np.diagonal(C, axis1=1, axis2=2), 0))
    return Result(gains, rgb.reshape(H, W, 3), x_reg, x_img, std.reshape(H, W, inv.N), skin_img, ing.reshape(H, W),
                  maha.reshape(H, W), rgb_model.reshape(H, W, 3), res.reshape(H, W, 3),
                  {"prior": prior_label, "mask_chi2": mask_chi2, "reg_sigma": reg_sigma, "reg_lambda": reg_lambda,
                   "reg_keep": list(reg_keep)})


def write_exr(path, arr):
    import OpenEXR

    arr = np.ascontiguousarray(arr, dtype=np.float16)
    channels = {"RGB": arr} if arr.ndim == 3 else {"Y": arr}
    header = {"compression": OpenEXR.ZIP_COMPRESSION, "type": OpenEXR.scanlineimage}
    with OpenEXR.File(header, channels) as f:
        f.write(path)


EXPORT_NAMES = {
    "melanin_fraction": "melanin",
    "eumelanin_ratio": "eumelanin_ratio",
    "blood_fraction": "hemoglobin",
    "oxygen_saturation": "oxygenation",
    "epidermis_thickness_cm": "epidermis_thickness_cm",
    "scattering_scale": "scattering_scale",
}


def export(result: Result, folder, prefix=""):
    """Write all maps as EXR (half float, ZIP). Non-skin texels keep their (unreliable) values; use the mask."""
    os.makedirs(folder, exist_ok=True)
    maps, rel = result.maps, result.relative_sigma()
    written = []
    for k, n in enumerate(inv.NAMES):
        for suffix, arr in (("", maps[..., k]), ("_sigma", rel[..., k])):
            path = os.path.join(folder, f"{prefix}{EXPORT_NAMES[n]}{suffix}.exr")
            write_exr(path, arr)
            written.append(path)
    extra = {
        "skin_mask": result.skin.astype(float),
        "albedo_input_calibrated": result.albedo,
        "albedo_from_maps": result.rgb_model,
        "colour_residual_coeffs": result.residual,
    }
    for name, arr in extra.items():
        path = os.path.join(folder, f"{prefix}{name}.exr")
        write_exr(path, arr)
        written.append(path)
    return written
