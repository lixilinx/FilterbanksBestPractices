"""Periodic-sequence-modulated filterbank design. 
Rewritten with AI based on materials:
    1, https://ieeexplore.ieee.org/document/8304771.
    2, https://github.com/lixilinx/FilterbanksBestPractices/blob/main/FilterBankDesign.m. 
"""

from dataclasses import dataclass, replace

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import linalg, sparse

FloatArray = NDArray[np.float64]
SparseMatrix = sparse.coo_matrix

@dataclass
class FilterBank:
    """The filterbank object. 
    """
    T: int # modulation sequence period
    B: int # decimation ratio/block size/hop size
    tau0: int # system delay
    h: ArrayLike # analysis prototype filter
    g: ArrayLike # synthesis prototype filter
    Gamma: ArrayLike | sparse.spmatrix | None = None # modulation product matrix
    i: int | None = None # analysis modulation shift
    j: int | None = None # synthesis modulation shift
    w_cut: float | None = None # prototype cutoff in rad/sample 
    zeta: float | None = None # synthesis stop-band weight
    symmetry: ArrayLike | None = None # symmetry[0]=1: g=flip(h)/MIRROR; symmetry[0]=-1: g=h/SAME; symmetry[1]=1: h=flip(h); symmetry[2]=1: g=flip(g)
    momentum: int | None = None # stop-band frequency weighting: 0, 1, or 2

    def __post_init__(self) -> None:
        assert np.ndim(self.h) == 1 and np.ndim(self.g) == 1
        Lh, Lg = np.size(self.h), np.size(self.g)
        assert self.T > 0
        assert 0 < self.B <= self.T
        assert Lh >= self.B and Lg >= self.B, "h/g must have >= B taps"
        assert self.B - 1 <= self.tau0 <= Lh + Lg - self.B - 1, "infeasible delay" # tighter bounds than matlab code 
        if self.Gamma is not None:
            assert self.Gamma.shape == (self.T, self.T)
        if self.w_cut is not None:
            assert 0 <= self.w_cut < np.pi 
        if self.zeta is not None:
            assert self.zeta >= 0
        if self.symmetry is not None:
            symmetry = np.asarray(self.symmetry).reshape(-1)
            assert symmetry.size == 3
            assert symmetry[0] == 0 or Lh == Lg, "MIRROR/SAME symmetry requires equal filter length"
        if self.momentum is not None:
            assert self.momentum in (0, 1, 2), "only momentum 0, 1, and 2 implemented"


@dataclass(frozen=True)
class DesignResult:
    """Hold design result returned by func design_filter_bank. 
    """
    filter_bank: FilterBank
    cost: float
    reconstruction_error: float
    iterations: int
    converged: bool


def _normalize_filter_bank(filter_bank: FilterBank) -> FilterBank:
    """Fill defaults. 
    """
    T, B, tau0 = filter_bank.T, filter_bank.B, filter_bank.tau0
    h = np.asarray(filter_bank.h, dtype=float).reshape(-1).copy()
    g = np.asarray(filter_bank.g, dtype=float).reshape(-1).copy()
    Gamma = filter_bank.Gamma
    if Gamma is None:
        Gamma = sparse.eye(T, format="csr")
    elif sparse.issparse(Gamma):
        Gamma = Gamma.tocsr().astype(float, copy=True)
    else:
        Gamma = np.asarray(Gamma, dtype=float).copy()

    return replace(
        filter_bank,
        Gamma=Gamma,
        i=(-tau0) % T if filter_bank.i is None else filter_bank.i,
        j=0 if filter_bank.j is None else filter_bank.j,
        h=h, g=g,
        w_cut=np.pi/B if filter_bank.w_cut is None else filter_bank.w_cut, # w_cut=pi if B=1, bad! 
        zeta=1.0 if filter_bank.zeta is None else filter_bank.zeta,
        symmetry=np.zeros(3) if filter_bank.symmetry is None else np.asarray(filter_bank.symmetry, dtype=float).reshape(-1).copy(),
        momentum=0 if filter_bank.momentum is None else filter_bank.momentum,
    )


def _ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def _active_terms(t: int, tau: int, B: int, Lh: int, Lg: int) -> range:
    """Return frame indices whose filter taps contribute at (t, tau).
    """
    n0 = max(_ceil_div(t - tau, B), (t - Lg) // B + 1)
    n1 = min(t // B, _ceil_div(Lh + t  - tau, B) - 1)
    return range(n0, n1 + 1)


def _valid_t_tau(Gamma: FloatArray, T: int, B: int, Lh: int, Lg: int, shift_i: int, shift_j: int) -> list[tuple[int, int]]:
    """Return all the (t, tau) pairs where M(t, tau) is not a zero matrix. 
    The Matlab code is greatly simplify after re-ordering the loop indices. 
    """
    tau_count = Lh + Lg - 1
    gamma_cols = (-np.arange(Lh) - shift_i) % T
    row_hits = Gamma[:, gamma_cols] != 0.0
    valid: list[tuple[int, int]] = []
    for t in range(B):
        hit = np.zeros(tau_count, dtype=bool)
        for n in range((t - Lg) // B + 1, t // B + 1):
            first = t - n * B
            hit[first : first + Lh] |= row_hits[(first + shift_j) % T]
        valid.extend((t, int(tau)) for tau in np.flatnonzero(hit))
    return valid


def _constraint_masks(Gamma: FloatArray, T: int, B: int, Lh: int, Lg: int, shift_i: int, shift_j: int, 
                      valid_t_tau: list[tuple[int, int]], tau0: int) -> tuple[SparseMatrix, FloatArray]:
    """Return all the M/mask matrices and targets. Mask matrices are row-stacked. 
    """
    size = Lh + Lg
    mask_count = len(valid_t_tau) + 1 # extra mask matrix for term h'h - g'g 
    targets = np.zeros(mask_count, dtype=np.float64)
    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []

    for mask_index, (t, tau) in enumerate(valid_t_tau):
        offset = mask_index * size
        for n in _active_terms(t, tau, B, Lh, Lg):
            h_index = n * B + tau - t
            g_index = Lh + t - n * B
            gamma_row = (t - n * B + shift_j) % T
            gamma_col = (t - tau - n * B - shift_i) % T
            value = 0.5 * Gamma[gamma_row, gamma_col]
            if value != 0.0:
                rows.extend((offset + h_index, offset + g_index))
                cols.extend((g_index, h_index))
                values.extend((value, value))
        if tau == tau0:
            targets[mask_index] = 1.0

    # the last mask blkdiag(I, -I) for the energy-balance term 
    offset = len(valid_t_tau) * size
    rows.extend(range(offset, offset + size))
    cols.extend(range(size))
    values.extend([1.0] * Lh + [-1.0] * Lg)

    stacked = sparse.coo_matrix((values, (rows, cols)), shape=(mask_count * size, size), dtype=np.float64)
    return stacked, targets


def _stopband_energy_matrix(length: int, w_cut: float, momentum: int) -> FloatArray:
    """Return the Toeplitz stop-band energy matrix. 
    """
    first_row = np.empty(length, dtype=np.float64)
    lag = 1.0 - np.arange(2, length + 1, dtype=np.float64)
    if momentum == 0:
        first_row[0] = np.pi - w_cut
        first_row[1:] = -np.sin(lag * w_cut) / lag
    elif momentum == 1:
        first_row[0] = (np.pi**2 - w_cut**2) / 2.0
        first_row[1:] = -(np.cos(lag * w_cut) + w_cut * lag * np.sin(lag * w_cut) - np.power(-1.0, lag)) / lag**2
    else: # momentum=2
        first_row[0] = (np.pi**3 - w_cut**3) / 3.0
        first_row[1:] = -((w_cut**2 * lag**2 - 2.0) * np.sin(lag * w_cut) + 2.0 * w_cut * lag * np.cos(lag * w_cut) - 2.0 * np.pi * lag * np.power(-1.0, lag)) / lag**3
    return linalg.toeplitz(first_row)


def _fixed_hessian(filter_bank: FilterBank, eta: float, lambda_: float) -> FloatArray:
    """Return the fixed part of the Hessian. 
    """
    h, g = np.asarray(filter_bank.h), np.asarray(filter_bank.g)
    symmetry = np.asarray(filter_bank.symmetry)
    w_cut = float(filter_bank.w_cut)
    zeta = float(filter_bank.zeta)
    momentum = int(filter_bank.momentum)
    Lh, Lg = h.size, g.size
    size = Lh + Lg

    hessian = np.zeros((size, size), dtype=np.float64)
    hessian[:Lh, :Lh] = _stopband_energy_matrix(Lh, w_cut, momentum)
    hessian[Lh:, Lh:] = zeta * _stopband_energy_matrix(Lg, w_cut, momentum)

    h_tap, g_tap = np.arange(Lh), np.arange(Lg)
    h_mirror, g_mirror = Lh - 1 - h_tap, Lg - 1 - g_tap
    all_taps = np.arange(size)
    hessian[all_taps, all_taps] += lambda_

    if symmetry[0] != 0.0:
        partner = h_mirror if symmetry[0] > 0.0 else h_tap
        hessian[h_tap, h_tap] += eta
        hessian[Lh + h_tap, Lh + h_tap] += eta
        hessian[h_tap, Lh + partner] -= eta
        hessian[Lh + partner, h_tap] -= eta

    h_paired = h_tap[h_tap != h_mirror] # exclude central tap for odd length filter 
    g_paired = g_tap[g_tap != g_mirror]
    if symmetry[1] != 0.0:
        hessian[h_paired, h_paired] += eta
        hessian[h_paired, Lh - 1 - h_paired] -= eta
    if symmetry[2] != 0.0:
        hessian[Lh + g_paired, Lh + g_paired] += eta
        hessian[Lh + g_paired, Lh + Lg - 1 - g_paired] -= eta

    return hessian


def _evaluate(x: FloatArray, fixed_hessian: FloatArray, stacked_masks: SparseMatrix, targets: FloatArray, B: int, 
              eta: float) -> tuple[float, float, FloatArray, FloatArray]:
    """Evaluate a filterbank prototype filter design x = (h, g).
    """
    all_mx = (stacked_masks @ x).reshape(-1, x.size).T
    errors = x @ all_mx - targets
    reconstruction_error = float(errors[:-1] @ errors[:-1] / B) # exclude the last one (the energy-balance term) 
    cost = float(0.5 * x @ fixed_hessian @ x + 0.5 * eta * (errors @ errors))
    return cost, reconstruction_error, errors, all_mx


def _search_direction(x: FloatArray, fixed_hessian: FloatArray, errors: FloatArray, all_mx: FloatArray, eta: float) -> FloatArray:
    """Compute the Hessian-preconditioned search direction. 
    """
    gradient = fixed_hessian @ x + 2.0 * eta * (all_mx @ errors)
    system = all_mx @ all_mx.T
    system *= 4.0 * eta 
    system += fixed_hessian
    try:
        direction = linalg.solve(system, gradient, assume_a="pos", check_finite=False) # Cholesky 
    except linalg.LinAlgError as error:
        raise np.linalg.LinAlgError("Hessian not positive definite; increase lambda_") from error 
    return np.asarray(direction, dtype=np.float64)


def design_filter_bank(filter_bank: FilterBank, eta: float, 
                       lambda_: float = 0.0, max_iterations: int = 100, tolerance: float = 1e-5) -> DesignResult:
    """Design prototype filters. 
    eta: penalty coefficient for reconstruction and constraints. 
    lambda_: L2 regularization coefficient. 
    """
    fb = _normalize_filter_bank(filter_bank)
    h, g = np.asarray(fb.h), np.asarray(fb.g)
    Gamma = fb.Gamma
    T, B, tau0, shift_i, shift_j = map(int, (fb.T, fb.B, fb.tau0, fb.i, fb.j))
    Lh, Lg = h.size, g.size

    gamma_lookup = Gamma.toarray() if sparse.issparse(Gamma) else np.asarray(Gamma)

    valid_t_tau = _valid_t_tau(gamma_lookup, T, B, Lh, Lg, shift_i, shift_j)
    if sum(tau == tau0 for _, tau in valid_t_tau) < B:
        raise ValueError("the given (i, j)-shift pair is infeasible")

    stacked_masks, targets = _constraint_masks(gamma_lookup, T, B, Lh, Lg, shift_i, shift_j, valid_t_tau, tau0)
    fixed_hessian = _fixed_hessian(fb, eta, lambda_)
    x = np.concatenate((h, g))
    converged = False 
    iterations = 0
    step_candidates = (1/2, 1/4, 1/16, 1/256, 1/65536)

    for iteration in range(1, max_iterations + 1):
        cost, _, errors, all_mx = _evaluate(x, fixed_hessian, stacked_masks, targets, B, eta)
        direction = _search_direction(x, fixed_hessian, errors, all_mx, eta)

        candidate = x
        for step in step_candidates:
            candidate = x - step * direction
            candidate_cost, *_ = _evaluate(candidate, fixed_hessian, stacked_masks, targets, B, eta)
            if candidate_cost < cost:
                break

        x = candidate
        iterations = iteration
        if np.max(np.abs(direction)) < tolerance:
            converged = True
            break

    final_cost, final_reconstruction_error, *_ = _evaluate(x, fixed_hessian, stacked_masks, targets, B, eta)
    designed = replace(fb, h=x[:Lh].copy(), g=x[Lh:].copy())
    return DesignResult(filter_bank=designed, cost=final_cost, reconstruction_error=final_reconstruction_error,
                        iterations=iterations, converged=converged)
