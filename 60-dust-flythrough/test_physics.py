"""Unit tests on synthetic inputs with known answers (no downloads needed).

    python -m pytest -q test_physics.py      # or: python test_physics.py
"""

from __future__ import annotations

import numpy as np
from astropy.io import fits

import config as C


def _synthetic_map(tmp_path, rho: float = 1e-3, inner: float = 0.05, nside: int = 4):
    """A HEALPix-times-distance file shaped like the real one, constant density everywhere."""
    npix = 12 * nside**2
    bounds = np.geomspace(50.0, 500.0, 21).astype(np.float32)
    radii = np.sqrt(bounds[:-1] * bounds[1:]).astype(np.float32)
    mean = np.full((radii.size, npix), rho, dtype=np.float32)
    hdr = {"NSIDE": nside, "ORDERING": "NEST", "CUNIT": "E of Zhang, Green, and Rix (2023)"}
    h_mean = fits.ImageHDU(mean, name="MEAN")
    h_std = fits.ImageHDU(0.5 * mean, name="STD.")
    for h in (h_mean, h_std):
        for k, v in hdr.items():
            h.header[k] = v
    h_in = fits.ImageHDU(np.full(npix, inner, dtype=np.float32), name="MEAN OF INTEGRATED INNER 50.0 PC")
    h_in.header["NSIDE"], h_in.header["ORDERING"] = nside, "NEST"
    h_in.header["CTYPE"] = "inner density integrated within 50.0 pc"
    t_c = fits.BinTableHDU.from_columns([fits.Column(name="radial pixel centers", format="E", array=radii)])
    t_b = fits.BinTableHDU.from_columns([fits.Column(name="radial pixel boundaries", format="E", array=bounds)])
    path = tmp_path / "synthetic.fits"
    fits.HDUList([fits.PrimaryHDU(), h_mean, h_std, t_c, t_b, h_in]).writeto(path)
    return path, bounds


def test_constant_density_integrates_exactly(tmp_path):
    from dustmap import DustSphere
    path, bounds = _synthetic_map(tmp_path)
    s = DustSphere(path=path)
    d = np.array([60.0, 100.0, 237.3, 499.0])
    l = np.array([10.0, 123.0, 250.0, 300.0])
    b = np.array([-40.0, 0.0, 12.0, 85.0])
    expected = 0.05 + 1e-3 * (d - bounds[0])
    assert np.allclose(s.extinction(l, b, d), expected, rtol=1e-5)
    # inside the inner sphere the inner integral is scaled linearly
    assert np.isclose(s.extinction(10.0, 0.0, 25.0)[()], 0.05 * 25.0 / bounds[0], rtol=1e-5)
    # beyond the last boundary the value saturates at the full column
    assert np.isclose(s.extinction(10.0, 0.0, 900.0)[()], 0.05 + 1e-3 * (bounds[-1] - bounds[0]), rtol=1e-5)
    # a constant field interpolates to itself in log space
    assert np.allclose(s.density_lbd(l, b, np.array([70.0, 150.0, 300.0, 420.0])), 1e-3, rtol=1e-5)


def test_coordinates_round_trip():
    from dustmap import lbd_to_xyz, xyz_to_lbd
    rng = np.random.default_rng(1)
    l, b, d = rng.uniform(0, 360, 1000), rng.uniform(-89, 89, 1000), rng.uniform(1, 2000, 1000)
    l2, b2, d2 = xyz_to_lbd(*lbd_to_xyz(l, b, d))
    assert np.allclose(l2, l) and np.allclose(b2, b) and np.allclose(d2, d)
    x, y, z = lbd_to_xyz(0.0, 0.0, 100.0)               # l = 0 points to the Galactic Centre (+x)
    assert np.isclose(x, 100.0) and abs(y) < 1e-9 and abs(z) < 1e-9


def test_unit_conversions():
    assert np.isclose(C.ZGR23_TO_NH, 1653.0, rtol=2e-3)  # O'Neill et al. (2024), eq. 4
    assert C.ZGR23_TO_AV == 2.8


def test_konietzka_wavelength_is_inverse_phase_derivative():
    import wave
    th = wave.KON_PUBLISHED["clouds"]
    s = np.linspace(0, 2000, 4001)
    lam_numeric = 1.0 / np.gradient(s / (th["p"] - th["gamma"] * s), s)
    assert np.allclose(wave.kon_wavelength(th, s)[1:-1], lam_numeric[1:-1], rtol=1e-4)
    # their quoted numbers: ~4 kpc near CMa OB1 (s = 0), ~0.4 kpc near Cygnus X (s ~ 2.27 kpc)
    assert np.isclose(wave.kon_wavelength(th, 0.0), th["p"])
    assert 380 < wave.kon_wavelength(th, 2274.0) < 420


def test_quadratic_baseline_passes_through_anchors():
    import wave
    th = wave.KON_PUBLISHED["clouds"]
    xa, ya = (th["x0"], th["x1"], th["x2"]), (th["y0"], th["y1"], th["y2"])
    y, _ = wave._quad(xa, ya, np.array(xa))
    assert np.allclose(y, ya)


def test_bubble_peak_finder_on_a_synthetic_profile():
    import bubble
    d = bubble.D_GRID
    # a wall: Gaussian bump at 150 pc (sigma 10 pc) on a faint floor, a second cloud at 300 pc
    prof = 1e-6 + 2e-4 * np.exp(-0.5 * ((d - 150) / 10) ** 2) + 5e-4 * np.exp(-0.5 * ((d - 300) / 8) ** 2)
    from scipy.ndimage import gaussian_filter1d
    out = bubble._first_peak(gaussian_filter1d(prof, bubble.SIGMA_SMOOTH), prof)
    d_peak, d_in, d_out = out[0], out[1], out[2]
    assert abs(d_peak - 150) <= 1
    # half-prominence width of the smoothed bump: FWHM of a sigma = sqrt(10^2 + 7^2) Gaussian
    fwhm = 2.3548 * np.hypot(10, 7)
    assert abs((d_out - d_in) - fwhm) < 2.0
    assert out[-1] == 2                                  # both clouds are found


def test_superbubble_radius_at_zucker_parameters():
    import history
    r = history.radius(14.39, 2.71, 1.06)                 # their medians: should give ~165 pc today
    assert 160 < r < 170


def test_viewer_quantisation_round_trip():
    import viewer
    rho = np.geomspace(3e-5, 0.15, 500)
    code = viewer.quantise(rho)
    back = 10 ** (viewer.LO + (code.astype(float) - 1) / 254 * (viewer.HI - viewer.LO))
    step = (viewer.HI - viewer.LO) / 254
    assert np.all(np.abs(np.log10(back / rho)) <= step / 2 + 1e-9)
    assert viewer.quantise(np.array([1e-6, 0.0]))[0] == 0


if __name__ == "__main__":
    import pathlib
    import tempfile
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            if "tmp_path" in fn.__code__.co_varnames[: fn.__code__.co_argcount]:
                with tempfile.TemporaryDirectory() as tmp:
                    fn(pathlib.Path(tmp))
            else:
                fn()
            print("ok", name)
