"""Smoke tests for the submission package (under a minute on one core).

    python3 -m pytest code/tests -q        (from ncs-submission/)
    python3 code/tests/test_smoke.py       (without pytest)

1. every solver module imports;
2. stored profiles solve the marching equations to the residual of the record and are smooth (m = 2);
3. the exact sparse Jacobian of the independent global solver matches central differences;
4. the Supplementary Tables are regenerated from data/outputs byte for byte;
5. one main and one supplementary figure are regenerated from data/outputs."""
import os
import subprocess
import sys
import tempfile
import filecmp

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.abspath(os.path.join(HERE, "..", ".."))
LIB = os.path.join(PKG, "code", "lib")
PROFILES = os.path.join(PKG, "data", "profiles")
sys.path.insert(0, LIB)


def test_imports():
    import importlib
    for mod in ("bq_logpolar", "bq_solver", "bq_newton", "bq_regrid", "bq_local_eig", "bq_stability", "bq_global",
                "ipm_solver", "ipm_local_eig", "ipm_wkb", "ipm_stability", "ipm_global"):
        importlib.import_module(mod)


def test_stored_profiles():
    """λ1 and λ3: the Newton residual stays at the stored level and m − 2 at the round-off of the 10-digit λ."""
    from ipm_solver import IPM
    from bq_newton import residual
    for n, lam in ((1, 0.4721297348), (3, 0.2415663353)):
        Y = np.load(os.path.join(PROFILES, f"ipm_profile_n{n}_lam{lam:.10f}_Nb32_hs0.0125.npy"))
        B = IPM(lam, hs=0.0125, Nb=32, s_sw=12.0, s_start=-20.0)
        R, info = residual(B, Y)
        assert R is not None
        assert np.abs(R).max() < 1e-9, (n, np.abs(R).max())
        assert abs(info["m"] - 2) < 1e-8, (n, info["m"] - 2)


def test_global_jacobian():
    from ipm_global import IPMGlobal
    G = IPMGlobal(s_min=-6.0, s_max=6.0, hs=0.1, Nb=12)
    rng = np.random.default_rng(0)
    U = np.concatenate([1e-2 * rng.standard_normal(3 * G.N), [0.3]])
    v = rng.standard_normal(U.size)
    eps = 1e-6
    fd = (G.residual(U + eps * v) - G.residual(U - eps * v)) / (2 * eps)
    err = np.abs(G.jacobian(U) @ v - fd).max() / np.abs(fd).max()
    assert err < 1e-8, err


def test_tables_reproduce():
    ref = os.path.join(PKG, "supplementary", "tables")
    with tempfile.TemporaryDirectory() as tmp:
        env = dict(os.environ, NCS_TABLES=tmp, MPLBACKEND="Agg")
        subprocess.run([sys.executable, os.path.join(PKG, "code", "analysis", "make_si_tables.py")], env=env,
                       check=True, capture_output=True)
        names = sorted(f for f in os.listdir(ref) if f.endswith(".tex"))
        assert names == sorted(os.listdir(tmp))
        match, mismatch, errors = filecmp.cmpfiles(ref, tmp, names, shallow=False)
        assert not mismatch and not errors, (mismatch, errors)


def test_figures_run():
    with tempfile.TemporaryDirectory() as tmp:
        env = dict(os.environ, NCS_FIGS=tmp, MPLBACKEND="Agg")
        subprocess.run([sys.executable, os.path.join(PKG, "code", "figures", "make_figures.py"), "fig2", "sfig1"],
                       env=env, check=True, capture_output=True)
        for f in ("Fig2.pdf", "Fig2.png", "SupplementaryFig1.pdf", "SupplementaryFig1.png"):
            assert os.path.getsize(os.path.join(tmp, f)) > 10_000, f


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
