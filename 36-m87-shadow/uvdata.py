"""Read EHT uvfits files and turn visibilities into the quantities we fit.

An interferometer does not record an image. Each pair of telescopes (a baseline)
records one complex number, the *visibility*, every few seconds: one Fourier
component of the sky brightness, at a spatial frequency (u, v) set by the
telescopes' separation as seen from the source, in wavelengths.

Conventions (the same as ehtim, the EHT's own library, verified in the tests):

    V(u, v) = integral I(x, y) exp(+2 pi i (u x + v y)) dx dy

with x the offset towards the east (increasing right ascension) and y towards
the north; position angles are measured east of north.

Station-based errors corrupt each visibility as V_ij -> g_i g_j* V_ij. Two
combinations cancel the unknown gains g exactly:

* closure phase, around a triangle:  arg(V_ij V_jk V_ki)
* log closure amplitude, around a quadrangle:  ln |V_ij| + ln |V_kl| - ln |V_ik| - ln |V_jl|

We build a *minimal* (non-redundant) set of each per scan, preferring the
highest-S/N combinations, and keep their full covariance matrix, as M87 Paper VI
does (section 2 and 4.1).
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
from astropy.io import fits

from config import INTRA_SITE_PAIRS, INTRA_SITE_UVDIST, MIN_SNR, SYSTEMATIC_FRACTION

SCAN_GAP_HOURS = 0.0165  # ehtim's default: a gap longer than ~59 s starts a new scan


@dataclass
class VisData:
    """Stokes I visibilities of one observation (one source, day and band)."""

    label: str
    source: str
    freq_hz: float
    mjd: int
    time: np.ndarray  # hours since 0h UT on day `mjd`
    t1: np.ndarray  # station codes
    t2: np.ndarray
    u: np.ndarray  # wavelengths
    v: np.ndarray
    vis: np.ndarray  # complex, Jy
    sigma: np.ndarray  # thermal noise of the complex visibility, Jy (per component)
    scan: np.ndarray  # integer scan index
    station_xyz: dict = field(default_factory=dict)  # geocentric metres
    ra_deg: float = np.nan
    dec_deg: float = np.nan

    def __len__(self) -> int:
        return len(self.vis)

    @property
    def uvdist(self) -> np.ndarray:
        return np.hypot(self.u, self.v)

    @property
    def amp(self) -> np.ndarray:
        return np.abs(self.vis)

    @property
    def snr(self) -> np.ndarray:
        return np.abs(self.vis) / self.sigma

    @property
    def stations(self) -> list[str]:
        return sorted(set(self.t1) | set(self.t2))

    def subset(self, mask: np.ndarray) -> "VisData":
        return replace(
            self,
            time=self.time[mask], t1=self.t1[mask], t2=self.t2[mask], u=self.u[mask], v=self.v[mask],
            vis=self.vis[mask], sigma=self.sigma[mask], scan=self.scan[mask],
        )

    def intra_site(self) -> np.ndarray:
        """True for baselines between co-located telescopes (ALMA-APEX, SMA-JCMT)."""
        pairs = np.array([frozenset((a, b)) in INTRA_SITE_PAIRS for a, b in zip(self.t1, self.t2)])
        return pairs | (self.uvdist < INTRA_SITE_UVDIST)


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------
def read_uvfits(path: Path | str, label: str | None = None) -> VisData:
    """Read an EHT Stokes-I uvfits file.

    The EHT releases store Stokes I in both parallel hands, RR = LL = I, with each
    hand's error set to sqrt(2) sigma_I (release README). So I = (RR + LL)/2 and
    sigma_I = 0.5 sqrt(sigma_RR^2 + sigma_LL^2), exactly what ehtim does.
    """
    path = Path(path)
    with fits.open(path) as hdul:
        grp = hdul[0].data
        hdr = hdul[0].header
        an = hdul["AIPS AN"]
        names = [n.strip() for n in an.data["ANNAME"]]
        nosta = list(an.data["NOSTA"])
        xyz = {n: np.array(p, dtype=float) for n, p in zip(names, an.data["STABXYZ"])}
        freq = float(an.header["FREQ"])

        # parameters are addressed by index: the two DATE parameters share a name
        parnames = list(grp.parnames)
        def par(name: str, occurrence: int = 0) -> np.ndarray:
            idx = [i for i, p in enumerate(parnames) if p == name][occurrence]
            return np.asarray(grp.par(idx), dtype=float)

        u_name = next(p for p in parnames if p.startswith("UU"))
        v_name = next(p for p in parnames if p.startswith("VV"))
        u = par(u_name) * freq  # light-seconds -> wavelengths
        v = par(v_name) * freq
        jd = par("DATE", 0) + par("DATE", 1)
        bl = np.round(par("BASELINE")).astype(int)
        a1, a2 = bl // 256, bl % 256

        data = np.asarray(grp["DATA"], dtype=float)  # (n, 1, 1, nif, nchan, nstokes, 3)
        data = data.reshape(data.shape[0], -1, data.shape[-2], 3)  # (n, nif*nchan, nstokes, 3)
        if data.shape[1] != 1:
            raise ValueError(f"{path.name}: expected one frequency channel, found {data.shape[1]}")
        rr = data[:, 0, 0, 0] + 1j * data[:, 0, 0, 1]
        ll = data[:, 0, 1, 0] + 1j * data[:, 0, 1, 1]
        w_rr, w_ll = data[:, 0, 0, 2], data[:, 0, 1, 2]

        scans_nx = None
        if "AIPS NX" in hdul:
            nx = hdul["AIPS NX"].data
            scans_nx = np.column_stack([nx["TIME"] - 0.5 * nx["TIME INTERVAL"],
                                        nx["TIME"] + 0.5 * nx["TIME INTERVAL"]]) * 24.0
        obj = str(hdr.get("OBJECT", "")).strip()
        ra, dec = float(hdr.get("OBSRA", np.nan)), float(hdr.get("OBSDEC", np.nan))

    good = (w_rr > 0) & (w_ll > 0) & np.isfinite(w_rr) & np.isfinite(w_ll) & np.isfinite(rr) & np.isfinite(ll)
    station_of = dict(zip(nosta, names))
    t1 = np.array([station_of[a] for a in a1])
    t2 = np.array([station_of[a] for a in a2])
    vis = 0.5 * (rr + ll)
    with np.errstate(divide="ignore"):
        sigma = 0.5 * np.sqrt(1.0 / w_rr + 1.0 / w_ll)

    mjd = int(np.floor(np.min(jd[good]) - 2400000.5))
    time = (jd - 2400000.5 - mjd) * 24.0
    out = VisData(
        label=label or path.stem, source=obj, freq_hz=freq, mjd=mjd, time=time[good], t1=t1[good],
        t2=t2[good], u=u[good], v=v[good], vis=vis[good], sigma=sigma[good],
        scan=np.zeros(good.sum(), dtype=int), station_xyz=xyz, ra_deg=ra, dec_deg=dec,
    )
    out.scan = assign_scans(out.time, scans_nx)
    return out


def assign_scans(time: np.ndarray, scan_table: np.ndarray | None = None) -> np.ndarray:
    """Scan index per record: from the NX table if present, else from gaps in time."""
    order = np.argsort(time)
    if scan_table is not None and len(scan_table):
        idx = np.full(len(time), -1, dtype=int)
        for k, (t0, t1) in enumerate(scan_table):
            idx[(time >= t0 - 1e-4) & (time <= t1 + 1e-4)] = k
        if np.all(idx >= 0):
            _, idx = np.unique(idx, return_inverse=True)
            return idx
    ts = time[order]
    new = np.concatenate([[0], (np.diff(ts) > SCAN_GAP_HOURS).astype(int)])
    idx_sorted = np.cumsum(new)
    idx = np.empty(len(time), dtype=int)
    idx[order] = idx_sorted
    return idx


# ---------------------------------------------------------------------------
# Averaging and amplitudes
# ---------------------------------------------------------------------------
def canonical_orientation(d: VisData, order: list[str] | None = None) -> VisData:
    """Flip baselines so t1 precedes t2 in `order` (alphabetical by default).

    Reversing a baseline conjugates its visibility and negates (u, v); the data
    are unchanged physically.
    """
    order = order or sorted(set(d.t1) | set(d.t2))
    rank = {s: i for i, s in enumerate(order)}
    flip = np.array([rank[a] > rank[b] for a, b in zip(d.t1, d.t2)])
    out = replace(d, t1=d.t1.copy(), t2=d.t2.copy(), u=d.u.copy(), v=d.v.copy(), vis=d.vis.copy())
    out.t1[flip], out.t2[flip] = d.t2[flip], d.t1[flip]
    out.u[flip], out.v[flip] = -d.u[flip], -d.v[flip]
    out.vis[flip] = np.conj(d.vis[flip])
    return out


def scan_average(d: VisData) -> VisData:
    """Coherently average each baseline over each scan.

    Exactly ehtim's ``avg_coherent(0, scan_avg=True)``: plain mean of the complex
    visibilities and of (u, v), with sigma = sqrt(sum sigma_i^2) / N. The EHT-HOPS
    pipeline has already removed the atmosphere's phase wander within each scan
    ("ad hoc phasing"), so the averaging is coherent; Paper III shows the
    coherence loss is negligible.
    """
    d = canonical_orientation(d)
    keys = np.array([f"{s}|{a}|{b}" for s, a, b in zip(d.scan, d.t1, d.t2)])
    uniq, inv = np.unique(keys, return_inverse=True)
    n = np.bincount(inv)

    def mean(x):
        return np.bincount(inv, weights=x) / n

    vis = mean(d.vis.real) + 1j * mean(d.vis.imag)
    sigma = np.sqrt(np.bincount(inv, weights=d.sigma**2)) / n
    first = np.array([np.flatnonzero(inv == k)[0] for k in range(len(uniq))])
    out = replace(
        d, time=mean(d.time), t1=d.t1[first], t2=d.t2[first], u=mean(d.u), v=mean(d.v), vis=vis,
        sigma=sigma, scan=d.scan[first],
    )
    order = np.lexsort((out.t2, out.t1, out.time))
    return out.subset(order)


def debiased_amplitude(amp: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    """Remove the upward noise bias of |V| (Rice distribution): sqrt(|V|^2 - sigma^2), floored at 0."""
    return np.sqrt(np.clip(amp**2 - sigma**2, 0.0, None))


def with_systematics(d: VisData, fraction: float = SYSTEMATIC_FRACTION) -> VisData:
    """Add a fractional 'non-closing' error in quadrature (1% in Paper VI)."""
    return replace(d, sigma=np.sqrt(d.sigma**2 + (fraction * np.abs(d.vis)) ** 2))


# ---------------------------------------------------------------------------
# Closure quantities
# ---------------------------------------------------------------------------
@dataclass
class ClosureSet:
    """A minimal set of closure quantities, as a linear map from per-baseline data.

    ``design`` (n_closures x n_vis) has +-1 entries: a closure phase is
    design @ (visibility phases), a log closure amplitude design @ ln|V|.
    ``cov`` is the full covariance of the observed closures (block-diagonal by scan).
    """

    kind: str  # "cphase" or "lcamp"
    design: np.ndarray
    values: np.ndarray  # radians or natural-log units
    cov: np.ndarray
    scan: np.ndarray
    time: np.ndarray
    stations: list  # tuple of station codes per closure
    _chol_inv: np.ndarray | None = None

    def __len__(self) -> int:
        return len(self.values)

    @property
    def sigma(self) -> np.ndarray:
        return np.sqrt(np.diag(self.cov))

    def whitener(self) -> np.ndarray:
        """Matrix W with W.T @ W = cov^-1, so chi^2 = |W r|^2."""
        if self._chol_inv is None:
            L = np.linalg.cholesky(self.cov)
            self._chol_inv = np.linalg.inv(L)
        return self._chol_inv

    def residual(self, model_values: np.ndarray) -> np.ndarray:
        r = self.values - model_values
        if self.kind == "cphase":
            r = np.angle(np.exp(1j * r))  # wrap into (-pi, pi]
        return r

    def chi2(self, model_values: np.ndarray) -> float:
        w = self.whitener() @ self.residual(model_values)
        return float(w @ w)


def _baseline_index(d: VisData) -> dict:
    return {(s, a, b): i for i, (s, a, b) in enumerate(zip(d.scan, d.t1, d.t2))}


def _row(n: int, entries) -> np.ndarray:
    r = np.zeros(n)
    for i, s in entries:
        r[i] += s
    return r


def _greedy_basis(rows: list[np.ndarray], sigmas: list[float], tol: float = 1e-9) -> list[int]:
    """Indices of a maximal linearly independent subset, taking the most precise first.

    This mirrors Paper VI's prescription (start from the redundant set and keep the
    highest-S/N closures up to the rank of the full set).
    """
    chosen: list[int] = []
    basis: list[np.ndarray] = []  # orthonormal rows (Gram-Schmidt)
    for k in np.argsort(sigmas, kind="stable"):
        r = rows[k].astype(float).copy()
        for b in basis:
            r -= (r @ b) * b
        nrm = np.linalg.norm(r)
        if nrm > tol:
            basis.append(r / nrm)
            chosen.append(int(k))
    return chosen


def closure_phases(d: VisData, *, exclude_intra_site: bool = True, min_snr: float = MIN_SNR) -> ClosureSet:
    """Minimal set of closure phases per scan, with full covariance.

    Triangles containing an intra-site baseline are dropped: their closure phase is
    trivially ~0 because the two co-located stations see the same Fourier point
    (Paper VI section 2.2). The phase error of each visibility is sigma/|V| (radians).
    """
    d = canonical_orientation(d)
    use = d.snr >= min_snr
    if exclude_intra_site:
        use &= ~d.intra_site()
    idx = _baseline_index(d)
    n = len(d)
    phase = np.angle(d.vis)
    var_phase = (d.sigma / np.abs(d.vis)) ** 2
    rows_all, vals, scans, times, sts = [], [], [], [], []
    for s in np.unique(d.scan):
        present = {(a, b) for (ss, a, b), i in idx.items() if ss == s and use[i]}
        st = sorted({x for ab in present for x in ab})
        rows, sig, trip = [], [], []
        for a, b, c in itertools.combinations(st, 3):
            if (a, b) in present and (b, c) in present and (a, c) in present:
                # V_ab V_bc V_ca = V_ab V_bc conj(V_ac)
                ent = [(idx[(s, a, b)], 1), (idx[(s, b, c)], 1), (idx[(s, a, c)], -1)]
                rows.append(_row(n, ent))
                sig.append(np.sqrt(sum(var_phase[i] for i, _ in ent)))
                trip.append((a, b, c))
        for k in _greedy_basis(rows, sig):
            rows_all.append(rows[k])
            scans.append(s)
            times.append(d.time[np.flatnonzero(d.scan == s)[0]])
            sts.append(trip[k])
    D = np.array(rows_all) if rows_all else np.zeros((0, n))
    values = np.angle(np.exp(1j * (D @ phase))) if len(D) else np.zeros(0)
    cov = D @ np.diag(var_phase) @ D.T if len(D) else np.zeros((0, 0))
    return ClosureSet("cphase", D, values, cov, np.array(scans), np.array(times), sts)


def log_closure_amplitudes(d: VisData, *, exclude_intra_site: bool = True,
                           min_snr: float = MIN_SNR) -> ClosureSet:
    """Minimal set of log closure amplitudes per scan, with full covariance.

    Uses noise-debiased amplitudes; the log-amplitude error of each baseline is
    sigma/A (Paper VI eq. 11-14). Amplitudes below `min_snr` are not used, since the
    Gaussian approximation for ln A fails at low S/N.
    """
    d = canonical_orientation(d)
    amp = debiased_amplitude(np.abs(d.vis), d.sigma)
    use = (amp > 0) & (amp / d.sigma >= min_snr)
    if exclude_intra_site:
        use &= ~d.intra_site()
    idx = _baseline_index(d)
    n = len(d)
    with np.errstate(divide="ignore"):
        lnamp = np.where(amp > 0, np.log(np.where(amp > 0, amp, 1.0)), 0.0)
        var_ln = np.where(amp > 0, (d.sigma / np.where(amp > 0, amp, 1.0)) ** 2, np.inf)

    def key(s, a, b):
        return idx.get((s, a, b)) if (s, a, b) in idx else idx.get((s, b, a))

    rows_all, scans, times, sts = [], [], [], []
    for s in np.unique(d.scan):
        present = {frozenset((a, b)) for (ss, a, b), i in idx.items() if ss == s and use[i]}
        st = sorted({x for ab in present for x in ab})
        rows, sig, quads = [], [], []
        for quad in itertools.combinations(st, 4):
            i, j, k, l = quad
            # the three ways of splitting four stations into two pairs; any two are independent
            for (p1, p2), (q1, q2) in (
                (((i, j), (k, l)), ((i, k), (j, l))),
                (((i, j), (k, l)), ((i, l), (j, k))),
                (((i, k), (j, l)), ((i, l), (j, k))),
            ):
                need = [p1, p2, q1, q2]
                if all(frozenset(x) in present for x in need):
                    ent = [(key(s, *p1), 1), (key(s, *p2), 1), (key(s, *q1), -1), (key(s, *q2), -1)]
                    rows.append(_row(n, ent))
                    sig.append(np.sqrt(sum(var_ln[ii] for ii, _ in ent)))
                    quads.append((p1, p2, q1, q2))
        for kk in _greedy_basis(rows, sig):
            rows_all.append(rows[kk])
            scans.append(s)
            times.append(d.time[np.flatnonzero(d.scan == s)[0]])
            sts.append(quads[kk])
    D = np.array(rows_all) if rows_all else np.zeros((0, n))
    var_safe = np.where(np.isfinite(var_ln), var_ln, 0.0)
    values = D @ lnamp if len(D) else np.zeros(0)
    cov = D @ np.diag(var_safe) @ D.T if len(D) else np.zeros((0, 0))
    return ClosureSet("lcamp", D, values, cov, np.array(scans), np.array(times), sts)


def trivial_closure_phases(d: VisData) -> np.ndarray:
    """Closure phases of triangles containing an intra-site baseline (should be ~0).

    A systematic check of the calibration: any departure beyond the noise means
    non-closing errors (Paper III uses the same test).
    """
    d = canonical_orientation(d)
    idx = _baseline_index(d)
    out = []
    for s in np.unique(d.scan):
        present = {(a, b) for (ss, a, b) in idx if ss == s}
        st = sorted({x for ab in present for x in ab})
        for a, b, c in itertools.combinations(st, 3):
            if (a, b) in present and (b, c) in present and (a, c) in present:
                pairs = [frozenset((a, b)), frozenset((b, c)), frozenset((a, c))]
                if any(p in INTRA_SITE_PAIRS for p in pairs):
                    i, j, k = idx[(s, a, b)], idx[(s, b, c)], idx[(s, a, c)]
                    psi = np.angle(d.vis[i] * d.vis[j] * np.conj(d.vis[k]))
                    sig = np.sqrt(sum((d.sigma[q] / abs(d.vis[q])) ** 2 for q in (i, j, k)))
                    out.append((psi, sig))
    return np.array(out)


@dataclass
class Prepared:
    """Everything the fitters need from one observation."""

    raw: VisData  # 10-s data
    avg: VisData  # scan-averaged, systematics added
    cphase: ClosureSet
    lcamp: ClosureSet

    @property
    def label(self) -> str:
        return self.raw.label


def prepare(path: Path | str, label: str | None = None, *, systematic: float = SYSTEMATIC_FRACTION,
            min_snr: float = MIN_SNR, exclude_intra_site: bool = False) -> Prepared:
    """Read, scan-average, add the 1% systematic and form the closure sets.

    As in Paper VI: closure phases leave out triangles with an intra-site baseline,
    while log closure amplitudes keep them (they tie the long baselines to the total
    flux on the shortest ones), so the models carry a large-scale component.
    """
    raw = read_uvfits(path, label)
    avg = with_systematics(scan_average(raw), systematic)
    cp = closure_phases(avg, exclude_intra_site=True, min_snr=min_snr)
    lca = log_closure_amplitudes(avg, exclude_intra_site=exclude_intra_site, min_snr=min_snr)
    return Prepared(raw, avg, cp, lca)


# ---------------------------------------------------------------------------
# Sgr A*: interstellar scattering and intrinsic variability
# ---------------------------------------------------------------------------
class Johnson2018Kernel:
    """Ensemble-average diffractive scattering kernel towards the Galactic Centre.

    The 'dipole' model of Johnson et al. (2018, ApJ 865, 104), in the approximate
    analytic form used by the EHT (Psaltis et al. 2018; ehtim's ScatteringModel,
    whose formulas and default parameters this mirrors line by line: ehtim's own
    class fails on numpy >= 2.5). The image is blurred by ~1.380 x 0.703 mas at 1 cm,
    scaling as wavelength^2, major axis at PA 81.9 deg: ~23 x 12 uas at 1.3 mm. The
    screen's 800-km inner scale makes the kernel fall off more slowly than a Gaussian
    on the longest baselines. The kernel is real and positive, so dividing by it
    ('deblurring') leaves closure phases untouched.
    """

    def __init__(self, alpha=1.38, d_obs_screen_cm=2.82 * 3.086e21, d_src_screen_cm=5.53 * 3.086e21,
                 theta_maj_mas=1.380, theta_min_mas=0.703, pa_deg=81.9, ref_wavelength_cm=1.0, r_in_cm=800e5):
        from scipy import integrate, special as sps
        from scipy.optimize import brentq

        self.alpha, self.r_in, self.lam_ref = alpha, r_in_cm, ref_wavelength_cm
        self.M = d_obs_screen_cm / d_src_screen_cm
        fwhm_fac = np.sqrt(2.0 * np.log(2.0)) / np.pi
        mas = np.pi / 180.0 / 3600.0 / 1000.0
        qbar = (2.0 / sps.gamma((2.0 - alpha) / 2.0)
                * (r_in_cm**2 * (1.0 + self.M) / (fwhm_fac * (ref_wavelength_cm / (2.0 * np.pi)) ** 2)) ** 2
                * ((theta_maj_mas**2 + theta_min_mas**2) * mas**2))
        c0 = (ref_wavelength_cm / (2.0 * np.pi)) ** 2 * qbar * sps.gamma(1.0 - alpha / 2.0) / (8.0 * np.pi**2 * r_in_cm**2)
        A = theta_maj_mas / theta_min_mas
        self.phi0 = np.radians(90.0 - pa_deg)
        self.Amaj = (r_in_cm * (1.0 + self.M) * theta_maj_mas * mas
                     / (fwhm_fac * (ref_wavelength_cm / (2.0 * np.pi)) * 2.0 * np.pi)) ** 2
        self.Amin = (r_in_cm * (1.0 + self.M) * theta_min_mas * mas
                     / (fwhm_fac * (ref_wavelength_cm / (2.0 * np.pi)) * 2.0 * np.pi)) ** 2
        # dipole angular distribution: solve for its anisotropy parameter k_zeta
        f = lambda k: (sps.hyp2f1((alpha + 2) / 2, 0.5, 2.0, -k) / sps.hyp2f1((alpha + 2) / 2, 1.5, 2.0, -k) - A**2)
        k = brentq(f, 1e-6, 1e3)
        pref = 1.0 / (2.0 * np.pi * sps.hyp2f1((alpha + 2) / 2, 0.5, 1.0, -k))
        P = lambda ph: pref * (1.0 + k * np.sin(ph - self.phi0) ** 2) ** (-(alpha + 2.0) / 2.0)
        int_maj = integrate.quad(lambda ph: abs(np.cos(self.phi0 - ph)) ** alpha * P(ph), 0, 2 * np.pi, limit=250)[0]
        int_min = integrate.quad(lambda ph: abs(np.sin(self.phi0 - ph)) ** alpha * P(ph), 0, 2 * np.pi, limit=250)[0]
        b_pref = c0 * 2.0 ** (2.0 - alpha) * np.sqrt(np.pi) / (alpha * sps.gamma((alpha + 1.0) / 2.0))
        self.Bmaj, self.Bmin = b_pref * int_maj, b_pref * int_min

    def _D(self, r, lam, A, B):
        a = self.alpha
        q = 2.0 * A / (a * B)
        return (lam / self.lam_ref) ** 2 * B * q ** (-a / (2.0 - a)) * (
            (1.0 + q ** (2.0 / (2.0 - a)) * (r / self.r_in) ** 2) ** (a / 2.0) - 1.0)

    def __call__(self, u, v, wavelength_cm):
        x = np.asarray(u) * wavelength_cm / (1.0 + self.M)
        y = np.asarray(v) * wavelength_cm / (1.0 + self.M)
        r, phi = np.hypot(x, y), np.arctan2(y, x)
        dmaj = self._D(r, wavelength_cm, self.Amaj, self.Bmaj)
        dmin = self._D(r, wavelength_cm, self.Amin, self.Bmin)
        dphi = 0.5 * (dmaj + dmin) + 0.5 * (dmaj - dmin) * np.cos(2.0 * (phi - self.phi0))
        return np.exp(-0.5 * dphi)


_KERNEL = None


def scattering_kernel(d: VisData) -> np.ndarray:
    """Johnson et al. (2018) kernel at each visibility's (u, v)."""
    global _KERNEL
    if _KERNEL is None:
        _KERNEL = Johnson2018Kernel()
    return _KERNEL(d.u, d.v, 2.99792458e10 / d.freq_hz)


def deblur(d: VisData) -> VisData:
    k = scattering_kernel(d)
    return replace(d, vis=d.vis / k, sigma=d.sigma / k)


def variability_sigma(uvdist: np.ndarray, *, a4: float, b: float, u0: float = 1.0e9, c: float = 2.0) -> np.ndarray:
    """Sgr A* Paper IV eq. 8: the broken power law of excess visibility scatter.

    sigma_var^2 = a^2 (|u|/u0)^c / (1 + (|u|/u0)^(b+c)), normalised so that
    sigma_var(4 Glambda) = a4 (in units of the total flux, as for light-curve
    normalised data). Paper IV measures a4 = 1.9%, b = 2.4 and u0 < 1.3 Glambda;
    c (the short-baseline index) is unconstrained and only matters below ~1 Glambda.
    """
    x = uvdist / u0
    x4 = 4.0e9 / u0
    a2 = a4**2 * (1.0 + x4 ** (b + c)) / x4**c
    return np.sqrt(a2 * x**c / (1.0 + x ** (b + c)))


def with_variability(d: VisData, *, a4: float, b: float, u0: float = 1.0e9, c: float = 2.0) -> VisData:
    sv = variability_sigma(d.uvdist, a4=a4, b=b, u0=u0, c=c)
    return replace(d, sigma=np.sqrt(d.sigma**2 + sv**2))


def prepare_sgra(paths: list, label: str, *, a4: float, b: float, systematic: float = SYSTEMATIC_FRACTION,
                 min_snr: float = MIN_SNR, u0: float = 1.0e9, c: float = 2.0) -> list[Prepared]:
    """Sgr A* for display: scan-average, add the 1% systematic and Paper IV's variability budget, deblur."""
    out = []
    for p in paths:
        raw = read_uvfits(p, label)
        avg = scan_average(raw)
        avg = with_systematics(avg, systematic)
        avg = with_variability(avg, a4=a4, b=b, u0=u0, c=c)  # Paper IV measured it on undeblurred data
        avg = deblur(avg)
        cp = closure_phases(avg, exclude_intra_site=True, min_snr=min_snr)
        lca = log_closure_amplitudes(avg, exclude_intra_site=False, min_snr=min_snr)
        out.append(Prepared(raw, avg, cp, lca))
    return out


def prepare_sgra_raw(paths: list, label: str, *, min_snr: float = MIN_SNR) -> list[Prepared]:
    """Sgr A* for the fitted-noise analysis: scan-averaged, thermal errors only, not deblurred.

    The minimal closure sets are chosen with the thermal errors; the noise model then
    sets their covariance at every likelihood call (fitting.NoiseModelLikelihood).
    """
    out = []
    for p in paths:
        raw = read_uvfits(p, label)
        avg = scan_average(raw)
        cp = closure_phases(avg, exclude_intra_site=True, min_snr=min_snr)
        lca = log_closure_amplitudes(avg, exclude_intra_site=False, min_snr=min_snr)
        out.append(Prepared(raw, avg, cp, lca))
    return out
