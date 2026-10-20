"""An offline 3D explorer of our basin of attraction, the flow lines, and the measured motions.

Writes ``figures/laniakea_3d.html``: a self-contained Plotly page (no network
needed) in supergalactic coordinates. Drag to rotate, scroll to zoom, click
legend entries to show or hide layers, hover over a group for its numbers.
"""

from __future__ import annotations

import numpy as np

import catalogs
import config
import density as dn
import plotstyle as ps

STEP = 4.0  # h^-1 Mpc, grid for the basin surfaces
HALF = 150.0


def _sg_grid():
    ax = np.arange(-HALF, HALF + STEP / 2, STEP)
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    return X, Y, Z, catalogs.sg_to_gal(pts)


def _indicator(labels, label, pts_gal, shape):
    ind = (labels == label).astype(np.float32)
    inside = np.all(np.abs(pts_gal) < config.RECON_HALF - 2, axis=1)
    out = np.zeros(len(pts_gal))
    out[inside] = dn.sample_field(ind, pts_gal[inside])
    r = np.linalg.norm(pts_gal, axis=1)
    out[r > config.RECON_RMAX - 3] = 0.0
    return out.reshape(shape)


def _surface(ind, color, name, opacity, visible=True):
    """A basin's boundary as a triangle mesh (marching cubes), far smaller than a volume trace."""
    import plotly.graph_objects as go
    from skimage import measure

    verts, faces, _, _ = measure.marching_cubes(ind, level=0.5, spacing=(STEP, STEP, STEP))
    verts = verts - HALF
    verts = np.round(verts, 1)
    return go.Mesh3d(x=verts[:, 0], y=verts[:, 1], z=verts[:, 2], i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
                     color=color, opacity=opacity, name=name, showlegend=True, hoverinfo="name", flatshading=False,
                     visible=True if visible else "legendonly",
                     lighting=dict(ambient=0.7, diffuse=0.5, specular=0.1))


def _streamlines(v_gal, seeds_sg, *, step=1.5, nmax=400):
    """Follow the flow (internal field) downhill from each seed; returns lists of SG points."""
    lines = []
    for s0 in seeds_sg:
        p = catalogs.sg_to_gal(np.asarray(s0, float))
        pts = [s0]
        for _ in range(nmax):
            if np.linalg.norm(p) > config.RECON_RMAX - 5:
                break
            v = dn.sample_field(v_gal, p[None, :])[0]
            sp = np.linalg.norm(v)
            if sp < 20:
                break
            p = p + step * v / sp
            pts.append(catalogs.gal_to_sg(p))
        if len(pts) > 8:
            lines.append(np.array(pts))
    return lines


def write_explorer(res, art) -> str:
    import plotly.graph_objects as go

    rec = art["rec"]
    beta = art["vc"]["beta_grid"]
    key = f"{config.BASIN_SMOOTH_MAIN:g}"
    bas = art["basins"][key]
    rep = bas["report"]
    X, Y, Z, pts = _sg_grid()
    traces = []
    ind = _indicator(bas["labels"], rep["ours_label"], pts, X.shape)
    traces.append(_surface(ind, ps.C_BASIN, "our basin of attraction", 0.28))
    names = {"Shapley": (-124, 75, -4), "Perseus-Pisces": (49, -10, -13), "Coma": (0, 71, 10)}
    for b in rep["largest"][:8]:
        if b["is_ours"]:
            continue
        att = np.array(b["attractor_sg"])
        nm = min(names, key=lambda n: np.linalg.norm(att - np.array(names[n])))
        if np.linalg.norm(att - np.array(names[nm])) > 30:
            continue
        ind_b = _indicator(bas["labels"], b["label"], pts, X.shape)
        traces.append(_surface(ind_b, ps.MUTED, f"{nm} basin", 0.18, visible=False))
    # Flow lines of the internal field, smoothed like the basins.
    _, cell = dn.grid_axes()
    extra = np.sqrt(config.BASIN_SMOOTH_MAIN**2 - config.RECON_SMOOTH**2)
    v = dn.velocity_field(dn.gaussian_smooth(dn.field_for(rec, beta), extra, cell), beta, cell)
    rng = np.random.default_rng(3)
    seeds = rng.uniform(-110, 110, size=(260, 3))
    seeds = seeds[np.linalg.norm(seeds, axis=1) < 115]
    lines = _streamlines(v, seeds)
    lx, ly, lz = [], [], []
    for ln in lines:
        ln = np.round(ln[::2], 1)
        lx += list(ln[:, 0]) + [None]
        ly += list(ln[:, 1]) + [None]
        lz += list(ln[:, 2]) + [None]
    traces.append(go.Scatter3d(x=lx, y=ly, z=lz, mode="lines", line=dict(color=ps.C_STREAM, width=2),
                               opacity=0.55, name="flow lines (2MRS gravity)", hoverinfo="skip"))
    # Measured motions.
    g = art["groups"]
    s = g[g["map_ok"] & (g["r_z"] < 100)].copy()
    s = s.sample(n=min(7000, len(s)), random_state=1)
    vv = np.clip(s["v_wf"].to_numpy(), -1000, 1000)
    cs = [[i / (len(ps.DIVERGING) - 1), c] for i, c in enumerate(ps.DIVERGING)]
    hover = [f"PGC {p}<br>cz (CMB) {cz:,.0f} km/s<br>distance {d:.0f} Mpc<br>peculiar velocity {u:+,.0f} km/s<br>{m}"
             for p, cz, d, u, m in zip(s["1PGC"], s["V3k"], s["d_mpc"], s["v_wf"], s["method"])]
    traces.append(go.Scatter3d(x=s["SGX_z"].round(1), y=s["SGY_z"].round(1), z=s["SGZ_z"].round(1), mode="markers",
                               marker=dict(size=2.2, color=vv, colorscale=cs, cmin=-1000, cmax=1000,
                                           colorbar=dict(title="peculiar<br>velocity<br>(km/s)", len=0.5, x=0.98)),
                               text=hover, hoverinfo="text", name="Cosmicflows-4 groups (< 100 h-1 Mpc)",
                               visible="legendonly"))
    lm = art["landmarks"]
    traces.append(go.Scatter3d(x=lm["SGX"], y=lm["SGY"], z=lm["SGZ"], mode="markers+text", text=lm["name"],
                               textposition="top center", marker=dict(size=4, color=ps.INK),
                               textfont=dict(size=11, color=ps.INK), name="clusters", hoverinfo="text"))
    att = rep["ours_attractor_sg"]
    traces.append(go.Scatter3d(x=[0, att[0]], y=[0, att[1]], z=[0, att[2]], mode="markers+text",
                               text=["Milky Way", "our attractor"], textposition="bottom center",
                               marker=dict(size=[6, 8], color=[ps.INK, ps.C_BASIN], symbol=["diamond", "x"]),
                               name="Milky Way and attractor", hoverinfo="text"))
    fig = go.Figure(traces)
    ax = dict(range=[-HALF, HALF], backgroundcolor=ps.SURFACE, gridcolor=ps.GRID, zerolinecolor=ps.AXIS)
    fig.update_layout(
        title=dict(text=f"Laniakea from 2MRS gravity: our basin of attraction ({key} h<sup>-1</sup> Mpc smoothing), "
                        "supergalactic coordinates", x=0.01, font=dict(size=15)),
        scene=dict(xaxis=dict(title="SGX (h<sup>-1</sup> Mpc)", **ax), yaxis=dict(title="SGY", **ax),
                   zaxis=dict(title="SGZ", **ax), aspectmode="cube",
                   camera=dict(eye=dict(x=1.35, y=-1.35, z=0.8))),
        paper_bgcolor=ps.SURFACE, font=dict(family="DejaVu Sans, sans-serif", color=ps.INK),
        legend=dict(x=0.01, y=0.95, bgcolor="rgba(252,252,251,0.8)"), margin=dict(l=0, r=0, t=50, b=0),
        annotations=[dict(text="Drag to rotate, scroll to zoom; click legend entries to show layers (other basins, "
                               "the measured groups). Hover a group for its numbers.", x=0.01, y=0.01, xref="paper",
                          yref="paper", showarrow=False, font=dict(size=11, color=ps.INK_SECONDARY), xanchor="left")],
    )
    out = config.FIGURE_DIR / "laniakea_3d.html"
    fig.write_html(out, include_plotlyjs=True, full_html=True, config={"displaylogo": False})
    return str(out)
