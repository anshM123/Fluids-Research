"""
qtt.py — quantized tensor-train (QTT / MPS) analysis of fields on 2^n grids.

Scale-interleaved ordering (Gourianov et al. 2022 convention): site l carries the l-th binary digit
of every coordinate (2^d-dimensional local index), coarse → fine. The bond after site l separates
scales coarser/finer than 2^{-l}. Ranks are computed by sequential truncated SVD (TT-SVD) with a
global relative Frobenius tolerance (eps split evenly over bonds, Oseledets 2011).
"""
import numpy as np


def to_scale_tensor(f):
    """f: array of shape (2^n,)*d → tensor of shape (2^d,)*n with scale-interleaved ordering."""
    d = f.ndim
    n = int(round(np.log2(f.shape[0])))
    assert all(s == 2**n for s in f.shape)
    # index i_k = sum_l b_{k,l} 2^{n-1-l}; reshape each axis into n bits (coarse first)
    t = f.reshape(*([2] * (n * d)))
    # current axis order: (b_{1,0..n-1}, b_{2,0..n-1}, ...) → want (b_{1,0}, b_{2,0}, ..., b_{1,1}, b_{2,1}, ...)
    perm = [k * n + l for l in range(n) for k in range(d)]
    t = np.transpose(t, perm)
    return t.reshape(*([2**d] * n)), n, d


def from_scale_tensor(t, n, d):
    t = t.reshape(*([2] * (n * d)))
    perm = [k * n + l for l in range(n) for k in range(d)]
    inv = np.argsort(perm)
    t = np.transpose(t, inv)
    return t.reshape(*([2**n] * d))


def tt_svd(f, eps=1e-3, chi_max=None, return_cores=False):
    """TT-SVD of scale-interleaved tensor. Returns list of bond ranks (length n−1) and relative error."""
    t, n, d = to_scale_tensor(np.asarray(f, dtype=float))
    q = 2**d
    norm = np.linalg.norm(t)
    delta = eps * norm / np.sqrt(max(n - 1, 1))
    cores, ranks = [], []
    r = 1
    C = t.reshape(q, -1)
    for l in range(n - 1):
        C = C.reshape(r * q, -1)
        U, s, Vt = np.linalg.svd(C, full_matrices=False)
        # smallest rank with tail norm <= delta
        tail = np.sqrt(np.cumsum(s[::-1] ** 2))[::-1]     # tail[k] = ||s[k:]||
        keep = int(np.sum(tail > delta))
        keep = max(keep, 1)
        if chi_max is not None:
            keep = min(keep, chi_max)
        cores.append(U[:, :keep].reshape(r, q, keep))
        C = s[:keep, None] * Vt[:keep]
        r = keep
        ranks.append(keep)
    cores.append(C.reshape(r, q, 1))
    if return_cores:
        return ranks, cores, (n, d)
    return ranks


def tt_full(cores, n, d):
    t = cores[0]
    for c in cores[1:]:
        t = np.tensordot(t, c, axes=([-1], [0]))
    t = t.reshape(*([2**d] * n))
    return from_scale_tensor(t, n, d)


def schmidt_spectra(f):
    """Singular values across every scale bond (exact, no truncation) — the 'interscale entanglement' spectrum."""
    t, n, d = to_scale_tensor(np.asarray(f, dtype=float))
    q = 2**d
    out = []
    for l in range(1, n):
        s = np.linalg.svd(t.reshape(q**l, -1), compute_uv=False)
        out.append(s)
    return out


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    x = np.linspace(0, 2 * np.pi, 256, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    f = np.sin(3 * X) * np.cos(5 * Y) + 0.3 * np.cos(X + 2 * Y)
    ranks, cores, (n, d) = tt_svd(f, eps=1e-10, return_cores=True)
    g = tt_full(cores, n, d)
    print("smooth field ranks:", ranks, "rel err", np.linalg.norm(f - g) / np.linalg.norm(f))
    ranks = tt_svd(rng.standard_normal((256, 256)), eps=1e-2)
    print("white noise ranks:", ranks)
