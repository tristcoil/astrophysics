"""ML extension: teach a small neural network to classify galaxies like the volunteers.

A convolutional neural network (CNN) learns image features -- edges, then
arcs and blobs, then arms and bulges -- by adjusting ~0.4 million weights
until its outputs match examples. Here the examples are GZ2 galaxies and the
target is *what the volunteers said*: for each of the first four questions
(smooth or features? edge-on? bar? spiral arms?) the network predicts the
fraction of volunteers giving each answer. It is trained by maximising the
probability of the actual vote counts (a multinomial likelihood), so a galaxy
with 50 votes counts for more than one with 20, and a question nobody was
asked contributes nothing -- the decision tree is handled for free. This is
the approach of Zoobot (Walmsley et al. 2022, 2023), in miniature.

Then it is tested on galaxies it never saw, against the volunteers and
against Zoobot's published predictions for the same galaxies (Galaxy Zoo
DESI). Zoobot is a much larger network trained on deeper DESI Legacy Surveys
images and on a newer volunteer project; it predicts its own volunteers'
answers, which on deeper images see more disks.

Run with ``python run.py --cnn`` (images are fetched on the first run).
"""

from __future__ import annotations

import json
import time

import numpy as np

import cnn_sample
import fetch_data as fd
from config import DATA_DIR, FIGURE_DIR, SEED

CNN_DIR = DATA_DIR / "cnn"
CNN_DIR.mkdir(parents=True, exist_ok=True)
IMAGES_PATH = CNN_DIR / "images_64.npz"
MODEL_PATH = CNN_DIR / "model.pt"
HISTORY_PATH = CNN_DIR / "history.json"
PRED_PATH = CNN_DIR / "test_predictions.npz"

#: Input size: the central half of each 424-pixel GZ2 image (4.2 times the
#: galaxy's 90%-light radius across), averaged down to 64 x 64 pixels.
PIX = 64
CROP = 0.25

#: The questions learned, and their answers in catalogue order.
QUESTIONS = {
    "smooth-or-features": ("smooth", "features", "artifact"),
    "edge-on": ("edgeon", "notedgeon"),
    "bar": ("bar", "nobar"),
    "spiral": ("spiral", "nospiral"),
}
EPOCHS = 20
BATCH = 128


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def _load_one(path):
    from PIL import Image

    img = Image.open(path).convert("RGB")
    n = img.size[0]
    k = int(round(n * CROP))
    img = img.crop((k, k, n - k, n - k)).resize((PIX, PIX), Image.Resampling.BOX)
    return np.asarray(img, np.uint8)


def has_archived_image(c: dict) -> np.ndarray:
    """Which galaxies have an image in the Zenodo archive (from cached files after the first run)."""
    index = fd.gz2_image_index(quiet=True)
    return np.isin(fd.gz2_asset_ids(c["dr7objid"]), index["asset_id"])


def prepare_images(c: dict, *, offline: bool = False, quiet: bool = False) -> dict:
    """The sample, its split, labels, and the images as a uint8 array (cached)."""
    s = cnn_sample.choose(c, has_archived_image(c))
    ids = c["dr7objid"][s["index"]]
    if not offline:
        order = np.argsort(s["split"] != "test", kind="stable")  # test images first
        fd.fetch_gz2_images(ids[order], quiet=quiet)
    if IMAGES_PATH.exists():
        with np.load(IMAGES_PATH) as f:
            if np.array_equal(f["dr7objid"], ids):
                return {**s, "dr7objid": ids, "images": f["images"], "ok": f["ok"]}
    images = np.zeros((len(ids), PIX, PIX, 3), np.uint8)
    ok = np.zeros(len(ids), bool)
    t0 = time.time()
    for k, o in enumerate(ids):
        p = fd.gz2_image_path(o)
        if p.exists():
            try:
                images[k] = _load_one(p)
                ok[k] = True
            except OSError:
                pass
    if ok.all():  # only cache a complete set
        np.savez(IMAGES_PATH, dr7objid=ids, images=images, ok=ok)
    if not quiet:
        print(f"  CNN images: {ok.sum():,d} of {len(ids):,d} prepared in {time.time() - t0:.0f}s")
    return {**s, "dr7objid": ids, "images": images, "ok": ok}


def vote_counts(c: dict, index) -> np.ndarray:
    """Votes per answer for the learned questions: shape (n, 9)."""
    cols = [a for answers in QUESTIONS.values() for a in answers]
    return np.column_stack([c["n_" + a][index] for a in cols]).astype(np.float32)


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

def build_model():
    """Four convolution stages (32, 64, 128, 128 filters), each halving the image.

    Kept small -- about 0.4 million weights, a third of the arithmetic of two
    convolutions per stage -- so it trains on a CPU in under an hour. Batch
    normalisation steadies training, and dropout before the last layer limits
    memorising the training set.
    """
    import torch.nn as nn

    def block(cin, cout):
        return nn.Sequential(nn.Conv2d(cin, cout, 3, padding=1, bias=False), nn.BatchNorm2d(cout),
                             nn.ReLU(inplace=True), nn.MaxPool2d(2))

    n_out = sum(len(a) for a in QUESTIONS.values())
    return nn.Sequential(block(3, 32), block(32, 64), block(64, 128), block(128, 128),
                         nn.Conv2d(128, 128, 3, padding=1, bias=False), nn.BatchNorm2d(128),
                         nn.ReLU(inplace=True), nn.AdaptiveAvgPool2d(1), nn.Flatten(),
                         nn.Dropout(0.3), nn.Linear(128, n_out))


def _split_logits(logits):
    out, k = [], 0
    for answers in QUESTIONS.values():
        out.append(logits[:, k:k + len(answers)])
        k += len(answers)
    return out


def multinomial_nll(logits, counts):
    """-log P(votes | predicted fractions), per vote, summed over the questions."""
    import torch

    loss, k, n = 0.0, 0, counts.sum()
    for part in _split_logits(logits):
        m = part.shape[1]
        loss = loss - (counts[:, k:k + m] * torch.log_softmax(part, dim=1)).sum()
        k += m
    return loss / n.clamp(min=1.0)


def predict_fractions(model, x, *, tta: bool = True, batch: int = 512) -> np.ndarray:
    """Predicted vote fractions, averaged over the 8 rotations/flips of each image."""
    import torch

    model.eval()
    outs = []
    with torch.no_grad():
        for lo in range(0, len(x), batch):
            xb = x[lo:lo + batch]
            views = [xb]
            if tta:
                views = [torch.rot90(v, k, dims=(2, 3)) for v in (xb, torch.flip(xb, dims=(3,)))
                         for k in range(4)]
            acc = 0
            for v in views:
                acc = acc + torch.cat([torch.softmax(p, 1) for p in _split_logits(model(v))], 1)
            outs.append((acc / len(views)).numpy())
    return np.concatenate(outs)


def _augment(xb):
    """Random rotation by 90-degree steps, flip, and a shift of up to 4 pixels.

    A galaxy has no preferred orientation on the sky, so these change nothing
    about its true shape: the network must learn shapes, not orientations.
    """
    import torch

    k = int(torch.randint(0, 4, (1,)))
    xb = torch.rot90(xb, k, dims=(2, 3))
    if torch.rand(1) < 0.5:
        xb = torch.flip(xb, dims=(3,))
    dx, dy = (int(v) for v in torch.randint(-4, 5, (2,)))
    return torch.roll(xb, shifts=(dy, dx), dims=(2, 3))


def train(data: dict, counts: np.ndarray, *, quiet: bool = False) -> dict:
    import torch

    torch.manual_seed(SEED)
    torch.set_num_threads(min(8, torch.get_num_threads()))
    x = torch.from_numpy(data["images"]).permute(0, 3, 1, 2).float() / 255.0
    y = torch.from_numpy(counts)
    tr = np.flatnonzero((data["split"] == "train") & data["ok"])
    va = np.flatnonzero((data["split"] == "val") & data["ok"])
    model = build_model()
    opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=1e-4)
    steps = EPOCHS * int(np.ceil(len(tr) / BATCH))
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=3e-3, total_steps=steps)
    rng = np.random.default_rng(SEED)
    history, best = [], (np.inf, None)
    for epoch in range(EPOCHS):
        model.train()
        t0, tot, n = time.time(), 0.0, 0
        perm = rng.permutation(tr)
        for lo in range(0, len(perm), BATCH):
            b = perm[lo:lo + BATCH]
            loss = multinomial_nll(model(_augment(x[b])), y[b])
            opt.zero_grad()
            loss.backward()
            opt.step()
            sched.step()
            tot += float(loss) * len(b)
            n += len(b)
        model.eval()
        with torch.no_grad():
            vl = float(multinomial_nll(model(x[va]), y[va]))
        history.append({"epoch": epoch + 1, "train": tot / n, "val": vl,
                        "seconds": time.time() - t0})
        if vl < best[0]:
            best = (vl, {k: v.clone() for k, v in model.state_dict().items()})
        if not quiet:
            print(f"    epoch {epoch + 1:2d}: train {tot / n:.4f}  val {vl:.4f}  "
                  f"({time.time() - t0:.0f}s)", flush=True)
    model.load_state_dict(best[1])
    torch.save(model.state_dict(), MODEL_PATH)
    HISTORY_PATH.write_text(json.dumps(history))
    return {"model": model, "history": history}


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def _metrics(pred, truth, n_votes) -> dict:
    """How close predicted fractions are to the volunteers', and the floor.

    The floor is the scatter the volunteers' own vote fractions would show if
    the same galaxy were shown to a fresh set of the same number of
    volunteers: binomial noise, sqrt(p (1 - p) / N). No predictor can beat
    it on average.
    """
    ok = np.isfinite(pred) & np.isfinite(truth) & (n_votes > 0)
    p, t, n = pred[ok], truth[ok], n_votes[ok]
    floor = np.sqrt(np.mean(t * (1 - t) / n))
    confident = (t > 0.8) | (t < 0.2)
    return {"n": int(ok.sum()), "rmse": float(np.sqrt(np.mean((p - t) ** 2))),
            "bias": float(np.mean(p - t)), "pearson_r": float(np.corrcoef(p, t)[0, 1]),
            "vote_noise_floor": float(floor),
            "n_confident": int(confident.sum()),
            "accuracy_confident": float(np.mean((p[confident] > 0.5) == (t[confident] > 0.5)))}


def bias_trend(model, data: dict, c: dict, d: dict, zoo_features=None) -> dict:
    """Do networks inherit 'distant galaxies look smoother'? (Fig. 02's panel B galaxies.)

    Mean 'features' fraction for galaxies of one luminosity and size, against
    redshift: the volunteers' raw votes, Hart et al.'s debiased ones, Zoobot's
    predictions (all GZ2 galaxies it matched), and our network's predictions
    for the held-out (validation and test) galaxies of that kind.
    """
    import torch

    sel = (d["good"] & (d["M"] >= -21.8) & (d["M"] < -21.4) & (d["R50_kpc"] >= 3)
           & (d["R50_kpc"] < 6))
    held = np.flatnonzero(((data["split"] == "test") | (data["split"] == "val")) & data["ok"])
    x = torch.from_numpy(data["images"][held]).permute(0, 3, 1, 2).float() / 255.0
    pf = predict_fractions(model, x)[:, 1]
    cnn_pf = np.full(len(c["dr7objid"]), np.nan)
    cnn_pf[data["index"][held]] = pf
    edges = np.array([0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.0925])
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = sel & (d["z"] >= lo) & (d["z"] < hi)
        row = {"z": float(np.median(d["z"][m])), "n": int(m.sum()),
               "raw": float(np.nanmean(c["wf_features"][m])),
               "hart16": float(np.nanmean(c["p_features"][m]))}
        if zoo_features is not None:
            row["zoobot"] = float(np.nanmean(zoo_features[m]))
            row["n_zoobot"] = int(np.isfinite(zoo_features[m]).sum())
        mc = m & np.isfinite(cnn_pf)
        row["n_cnn"] = int(mc.sum())
        row["cnn"] = float(np.mean(cnn_pf[mc])) if mc.sum() >= 15 else None
        row["raw_on_cnn"] = float(np.nanmean(c["wf_features"][mc])) if mc.sum() >= 15 else None
        rows.append(row)
    return {"bins": rows}


def evaluate(model, data: dict, c: dict) -> dict:
    import pandas as pd
    import torch

    te = np.flatnonzero((data["split"] == "test") & data["ok"])
    x = torch.from_numpy(data["images"][te]).permute(0, 3, 1, 2).float() / 255.0
    frac = predict_fractions(model, x)
    gi = data["index"][te]
    cols = [a for answers in QUESTIONS.values() for a in answers]
    pred = {a: frac[:, k] for k, a in enumerate(cols)}
    np.savez(PRED_PATH, dr7objid=c["dr7objid"][gi], **{f"pred_{a}": pred[a] for a in cols})

    def nq(answers):
        return sum(c["n_" + a][gi] for a in answers).astype(float)

    truth = {a: c["wf_" + a][gi].astype(float) for a in cols}
    out = {"n_test": int(len(te)), "cnn": {}, "zoobot": {}}
    out["cnn"]["features"] = _metrics(pred["features"], truth["features"], nq(QUESTIONS["smooth-or-features"]))
    out["cnn"]["edgeon"] = _metrics(pred["edgeon"], truth["edgeon"], nq(QUESTIONS["edge-on"]))
    out["cnn"]["bar"] = _metrics(pred["bar"], truth["bar"], nq(QUESTIONS["bar"]))
    out["cnn"]["spiral"] = _metrics(pred["spiral"], truth["spiral"], nq(QUESTIONS["spiral"]))

    # Zoobot (GZ DESI), for the same test galaxies where they match.
    if fd.GZDESI_SUBSET_PATH.exists():
        z = pd.read_parquet(fd.GZDESI_SUBSET_PATH).drop_duplicates("dr7objid").set_index("dr7objid")
        zz = z.reindex(c["dr7objid"][gi])
        zf = zz["smooth-or-featured_featured-or-disk_fraction"].to_numpy(float)
        ze = zz["disk-edge-on_yes_fraction"].to_numpy(float)
        zb = (zz["bar_strong_fraction"] + zz["bar_weak_fraction"]).to_numpy(float)
        zbs = zz["bar_strong_fraction"].to_numpy(float)
        zs = zz["has-spiral-arms_yes_fraction"].to_numpy(float)
        m = np.isfinite(zf)
        out["zoobot"]["n_matched"] = int(m.sum())
        out["zoobot"]["features"] = _metrics(zf, truth["features"], nq(QUESTIONS["smooth-or-features"]))
        out["zoobot"]["edgeon"] = _metrics(ze, truth["edgeon"], nq(QUESTIONS["edge-on"]))
        out["zoobot"]["bar_strong_or_weak"] = _metrics(zb, truth["bar"], nq(QUESTIONS["bar"]))
        out["zoobot"]["bar_strong"] = _metrics(zbs, truth["bar"], nq(QUESTIONS["bar"]))
        out["zoobot"]["spiral"] = _metrics(zs, truth["spiral"], nq(QUESTIONS["spiral"]))
        # The CNN on exactly the matched galaxies, for a like-for-like comparison.
        out["cnn_on_matched"] = {"features": _metrics(np.where(m, pred["features"], np.nan),
                                                      truth["features"],
                                                      nq(QUESTIONS["smooth-or-features"]))}
        out["_zoobot_features"] = zf
    out["_pred"] = pred
    out["_truth"] = truth
    out["_gi"] = gi
    return out


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def figure(ev: dict, history: list, c: dict) -> list:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    from plotstyle import (AQUA, BLUE, BLUE_RAMP, C_OURS, INK, INK_SECONDARY, MUTED, ORANGE, SURFACE,
                           caption, save, use_style)

    use_style()
    cmap = LinearSegmentedColormap.from_list("blues", BLUE_RAMP)
    fig = plt.figure(figsize=(13, 8.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.5, wspace=0.34)

    ax = fig.add_subplot(gs[0, 0])
    ep = [h["epoch"] for h in history]
    ax.plot(ep, [h["train"] for h in history], color=MUTED, lw=1.6, label="training galaxies")
    ax.plot(ep, [h["val"] for h in history], color=BLUE, lw=1.8, label="validation galaxies")
    ax.set_xlabel("epoch (one pass through the training set)")
    ax.set_ylabel("loss (−log probability of\nthe votes, per vote)")
    ax.legend(loc="upper right")
    ax.set_title("A. Learning", fontsize=10)

    truth, pred = ev["_truth"]["features"], ev["_pred"]["features"]
    for k, (yv, name) in enumerate(((pred, "our CNN (SDSS images)"),
                                     (ev.get("_zoobot_features"), "Zoobot (deeper DESI images)"))):
        ax = fig.add_subplot(gs[0, 1 + k])
        if yv is None:
            ax.axis("off")
            continue
        ok = np.isfinite(yv) & np.isfinite(truth)
        ax.hexbin(truth[ok], yv[ok], gridsize=28, extent=(0, 1, 0, 1), bins="log", mincnt=1,
                  cmap=cmap, linewidths=0, rasterized=True)
        ax.plot([0, 1], [0, 1], color=INK, lw=0.8)
        key = "cnn" if k == 0 else "zoobot"
        mt = ev[key]["features"]
        ax.text(0.03, 0.97, f"rms difference {mt['rmse']:.3f}\nmean offset {mt['bias']:+.3f}\n"
                            f"r = {mt['pearson_r']:.3f}", transform=ax.transAxes, va="top",
                fontsize=8, color=INK)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("GZ2 volunteers: fraction saying 'features'")
        ax.set_ylabel(f"{name.split(' (')[0]}: predicted fraction")
        ax.set_title(("B. " if k == 0 else "C. ") + name, fontsize=10)

    ax = fig.add_subplot(gs[1, 0])
    qs = [("features", "features\nor disk"), ("edgeon", "edge-on"), ("bar", "bar"),
          ("spiral", "spiral\narms")]
    xq = np.arange(len(qs))
    cnn_rmse = [ev["cnn"][q]["rmse"] for q, _ in qs]
    floor = [ev["cnn"][q]["vote_noise_floor"] for q, _ in qs]
    ax.bar(xq - 0.2, cnn_rmse, width=0.38, color=C_OURS, label="our CNN")
    ax.bar(xq + 0.2, floor, width=0.38, color=MUTED, label="volunteers' own noise")
    ax.set_xticks(xq)
    ax.set_xticklabels([l for _, l in qs], fontsize=8)
    ax.set_ylabel("rms difference from the\nvolunteers' vote fraction")
    ax.legend(loc="upper left")
    ax.set_title("D. Error against the noise floor", fontsize=10)
    ax.grid(axis="x", visible=False)

    ax = fig.add_subplot(gs[1, 1])
    acc = [100 * ev["cnn"][q]["accuracy_confident"] for q, _ in qs]
    zq = {"features": "features", "edgeon": "edgeon", "bar": "bar_strong_or_weak", "spiral": "spiral"}
    zacc = [100 * ev["zoobot"][zq[q]]["accuracy_confident"] if ev.get("zoobot", {}).get(zq[q]) else np.nan
            for q, _ in qs]
    ax.bar(xq - 0.2, acc, width=0.38, color=C_OURS, label="our CNN")
    ax.bar(xq + 0.2, zacc, width=0.38, color=ORANGE, label="Zoobot")
    for xx, a, b in zip(xq, acc, zacc):
        ax.text(xx - 0.2, a + 1, f"{a:.0f}", ha="center", fontsize=7.5)
        if np.isfinite(b):
            ax.text(xx + 0.2, b + 1, f"{b:.0f}", ha="center", fontsize=7.5)
    ax.set_xticks(xq)
    ax.set_xticklabels([l for _, l in qs], fontsize=8)
    ax.set_ylim(0, 105)
    ax.set_ylabel("% agreeing with a clear\nvolunteer majority (>80%)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    ax.set_title("E. Where the volunteers were sure", fontsize=10)
    ax.grid(axis="x", visible=False)

    ax = fig.add_subplot(gs[1, 2])
    bt = ev.get("bias_trend")
    if bt:
        rows = bt["bins"]
        z = np.array([r["z"] for r in rows])
        for key, col, lab in (("raw", AQUA, "GZ2 volunteers, raw"), ("zoobot", ORANGE, "Zoobot"),
                              ("hart16", BLUE, "Hart et al. debiased")):
            y = np.array([r.get(key, np.nan) for r in rows], float)
            ax.plot(z, y, color=col, marker="o", ms=4, lw=1.6, mec=SURFACE, mew=1, label=lab)
        zc = np.array([r["z"] for r in rows if r["cnn"] is not None])
        yc = np.array([r["cnn"] for r in rows if r["cnn"] is not None])
        yr = np.array([r["raw_on_cnn"] for r in rows if r["cnn"] is not None])
        if len(zc):
            ax.plot(zc, yc, color=INK, marker="D", ms=4.5, lw=0, label="our CNN, held-out galaxies")
            ax.plot(zc, yr, color=AQUA, marker="D", ms=4.5, lw=0, mfc="none", mew=1.2,
                    label="their raw votes")
        ax.set_ylim(0, 0.85)
        ax.set_xlabel("redshift z")
        ax.set_ylabel("mean 'features' fraction,\none kind of galaxy (Fig. 02B)")
        ax.legend(loc="lower left", fontsize=7.3)
    ax.set_title("F. Networks learn the bias too", fontsize=10)

    n_w = sum(p.numel() for p in build_model().parameters())
    caption(fig, f"A small convolutional network ({n_w / 1e6:.2f} million weights) trained on the votes "
                 "for 9,000 GZ2 galaxies (the central half of each volunteer image at 64 x 64 pixels, "
                 "randomly rotated and flipped) to predict the fraction of volunteers giving each answer "
                 "to the first four questions. Tested on ~2,000 galaxies it never saw. Zoobot: Walmsley "
                 "et al. (2023)'s published predictions from DESI Legacy Surveys images (its 'bar' = strong "
                 "+ weak, an answer GZ2 did not offer). D: the noise floor "
                 "is how much a vote fraction would change with a fresh set of the same number of "
                 "volunteers. F: galaxies of one luminosity and size at increasing distance; networks "
                 "trained to predict volunteers predict their biases too.", y=-0.01)
    return save(fig, FIGURE_DIR / "08_neural_network")


def _zoobot_features(c: dict):
    import pandas as pd

    if not fd.GZDESI_SUBSET_PATH.exists():
        return None
    z = pd.read_parquet(fd.GZDESI_SUBSET_PATH).drop_duplicates("dr7objid").set_index("dr7objid")
    return z.reindex(c["dr7objid"])["smooth-or-featured_featured-or-disk_fraction"].to_numpy(float)


def run(state: dict, *, offline: bool = False, quiet: bool = False, retrain: bool = False) -> dict:
    import torch

    c = state["c"]
    if not quiet:
        print("Neural network")
    data = prepare_images(c, offline=offline, quiet=quiet)
    counts = vote_counts(c, data["index"])
    if MODEL_PATH.exists() and HISTORY_PATH.exists() and not retrain:
        model = build_model()
        model.load_state_dict(torch.load(MODEL_PATH))
        history = json.loads(HISTORY_PATH.read_text())
    else:
        out = train(data, counts, quiet=quiet)
        model, history = out["model"], out["history"]
    if not offline:
        ra = np.where(np.isfinite(c["ra_sdss"]), c["ra_sdss"], c["ra"])
        dec = np.where(np.isfinite(c["dec_sdss"]), c["dec_sdss"], c["dec"])
        fd.fetch_gzdesi_subset(ra, dec, c["dr7objid"], quiet=quiet)
    ev = evaluate(model, data, c)
    ev["bias_trend"] = bias_trend(model, data, c, state["d"], _zoobot_features(c))
    figure(ev, history, c)
    state["cnn"] = ev
    n_params = sum(p.numel() for p in model.parameters())
    res = {k: v for k, v in ev.items() if not k.startswith("_")}
    res.update({"n_parameters": int(n_params), "n_train": int(((data["split"] == "train") & data["ok"]).sum()),
                "n_val": int(((data["split"] == "val") & data["ok"]).sum()),
                "epochs": len(history), "best_val_loss": float(min(h["val"] for h in history))})
    return res
