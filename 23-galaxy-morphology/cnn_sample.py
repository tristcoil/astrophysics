"""Which galaxies the neural network learns from, and which it is tested on.

Kept apart from ``cnn.py`` (and from PyTorch) so the images can be fetched
before any training code runs.

The split is fixed by a hash of each galaxy's ID, not by a random draw that a
code change could reshuffle: a galaxy is in the test set forever, and the
network never sees it during training or model selection.
"""

from __future__ import annotations

import hashlib

import numpy as np

#: How many galaxies, and how they split. Limited by how fast the image
#: archive serves files (~3 per second).
N_TOTAL = 12_000
SPLIT = {"train": 0.75, "val": 0.08, "test": 0.17}


def _bucket(objid: int) -> float:
    """A stable pseudo-random number in [0, 1) from the galaxy's ID."""
    h = hashlib.sha256(str(int(objid)).encode()).digest()
    return int.from_bytes(h[:8], "little") / 2**64


def choose(c: dict, has_image: np.ndarray, n_total: int = N_TOTAL) -> dict:
    """A random sample of well-classified GZ2 galaxies, split train/val/test.

    Every GZ2 galaxy with an archived image and at least 20 classifications is
    eligible; the ``n_total`` with the smallest hash values are taken, which is
    a uniform random sample that never changes.
    """
    ids = np.asarray(c["dr7objid"], np.int64)
    ok = has_image & (c["n_class"] >= 20)
    cand = np.flatnonzero(ok)
    keys = np.array([_bucket(o) for o in ids[cand]])
    chosen = cand[np.argsort(keys)[:n_total]]
    # A second, independent hash decides the split.
    split_key = np.array([_bucket(o * 7 + 3) for o in ids[chosen]])
    edges = np.cumsum([SPLIT["train"], SPLIT["val"]])
    part = np.where(split_key < edges[0], "train", np.where(split_key < edges[1], "val", "test"))
    return {"index": chosen, "split": part}
