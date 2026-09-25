"""Minimal numpy runner for the BioSkin (Aliaga et al. 2023) pretrained decoder.

Loads the PyTorch checkpoint without torch (zip + pickle of raw storages) and
evaluates the spectral decoder p -> R(lambda), 380-1000 nm / 2 nm (310 values).
Parameter warping and the RGB conversion follow the BioSkin repository
(MIT licence): https://github.com/facebookresearch/BioSkin
"""

import json
import pickle
import zipfile

import numpy as np

THICK_MIN, THICK_MAX = 0.001, 0.035
DTYPES = {"FloatStorage": np.float32, "DoubleStorage": np.float64, "LongStorage": np.int64}


class _Storage:
    def __init__(self, dtype, key):
        self.dtype, self.key = dtype, key


def load_state_dict(path):
    z = zipfile.ZipFile(path)
    root = z.namelist()[0].split("/")[0]

    class U(pickle.Unpickler):
        def find_class(self, module, name):
            if name in DTYPES:
                return name
            if name == "_rebuild_tensor_v2":
                def rebuild(storage, offset, size, stride, *args):
                    raw = np.frombuffer(z.read(f"{root}/data/{storage.key}"), dtype=storage.dtype)
                    return np.lib.stride_tricks.as_strided(raw[offset:], shape=size, strides=[s * raw.itemsize for s in stride]).copy()
                return rebuild
            if name == "OrderedDict":
                from collections import OrderedDict
                return OrderedDict
            return super().find_class(module, name)

        def persistent_load(self, pid):
            _, stype, key, _, _ = pid
            return _Storage(DTYPES[stype if isinstance(stype, str) else stype.__name__], key)

    return U(z.open(f"{root}/data.pkl")).load()


class BioSkinDecoder:
    def __init__(self, model_prefix):
        self.params = json.load(open(model_prefix + ".json"))
        self.w = {k.replace("module.", ""): v for k, v in load_state_dict(model_prefix + ".pt").items()}
        self.lam = np.arange(self.params["wavelength_begin"], self.params["wavelength_end"], self.params["spectral_resolution"])

    @staticmethod
    def unwarp(mel, blood, thick, so2, ratio):
        """Linear skin properties -> network input space (see bioskin/parameters/params_io.py)."""
        return np.stack([np.cbrt(mel), blood ** 0.25, (thick - THICK_MIN) / (THICK_MAX - THICK_MIN), so2, ratio], -1)

    def decode(self, x):
        w = self.w
        h = np.tanh(x @ w["fc_dec_in.weight"].T + w["fc_dec_in.bias"])
        h = np.tanh(h @ w["fc_dec.weight"].T + w["fc_dec.bias"])
        return 1 / (1 + np.exp(-(h @ w["fc_dec_out.weight"].T + w["fc_dec_out.bias"])))
