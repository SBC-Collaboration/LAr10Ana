from typing import NamedTuple

import numpy as np
from sklearn.cluster import DBSCAN, HDBSCAN


class Point(NamedTuple):
    """One bubble-finder detection. Every algo takes a list of these."""
    x: float             # pixel coordinates of the detection
    y: float
    frame: int           # frame number the detection came from
    significance: float  # finder's confidence, 0 to 1 (Hough vote ratio)
    radius: float        # radius estimate in pixels


def clusters_nearest(points, max_dist=15, min_frames=3):
    """Non-greedy: assign each point to its nearest cluster, not the first one in range."""
    # [cx, cy, set_of_frames]
    clusters = []
    for x, y, frame, *_ in sorted(points, key=lambda p: p.frame):
        # find the closest cluster center, and how far away it is
        best = None
        best_dist = max_dist**2
        for cluster in clusters:
            cx, cy, _ = cluster
            dist = (x - cx) ** 2 + (y - cy) ** 2
            if dist <= best_dist:
                best = cluster
                best_dist = dist

        if best is not None:
            # add this point's frame to the nearest cluster in range
            best[2].add(frame)
        else:
            # nothing close enough, start a new cluster
            clusters.append([x, y, {frame}])

    # clusters that appear in enough frames are considered bubbles
    return sum(1 for c in clusters if len(c[2]) >= min_frames)


def clusters_greedy(points, max_dist=15, min_frames=3):
    """Count clusters of points in the same (x, y) location across multiple
    frames. Accept first cluster(greedy) within max_dist(not necessarily the closest)"""
    # [cx, cy, set_of_frames]
    clusters = []
    for x, y, frame, *_ in sorted(points, key=lambda p: p.frame):
        # first cluster this point is close enough to, if any
        match = next(
            (c for c in clusters if (x - c[0]) ** 2 + (y - c[1]) ** 2 <= max_dist**2),
            None,
        )
        if match is not None:
            match[2].add(frame)
        else:
            # nothing close enough, start a new cluster
            clusters.append([x, y, {frame}])

    return sum(1 for c in clusters if len(c[2]) >= min_frames)


def clusters_track(points, max_dist=15, min_frames=3):
    """Track where the cluster was last seen, and then move it the new position"""
    # [cx, cy, set_of_frames]  cx, cy follow the most recent point in the cluster
    clusters = []
    # process in frame order so "last position" means the previous frame
    for x, y, frame, *_ in sorted(points, key=lambda p: p.frame):
        # find the closest cluster's current (last-seen) position
        best = None
        best_dist = max_dist**2
        for cluster in clusters:
            cx, cy, _ = cluster
            dist = (x - cx) ** 2 + (y - cy) ** 2
            if dist <= best_dist:
                best = cluster
                best_dist = dist

        if best is not None:
            # extend the track and move its center to where the bubble is now
            best[0], best[1] = x, y
            best[2].add(frame)
        else:
            # nothing close enough, start a new track
            clusters.append([x, y, {frame}])

    return sum(1 for c in clusters if len(c[2]) >= min_frames)


def clusters_significance(points, max_dist=15, min_frames=3, min_sig=0.5):
    """Drop low-confidence detections (significance < min_sig), then track-cluster the rest.
    Significance is the bubble finder's Hough vote ratio: faint/noisy hits score low."""
    strong = [p for p in points if p.significance >= min_sig]
    return clusters_track(strong, max_dist, min_frames)


def clusters_radius(points, max_dist=15, min_frames=3):
    """Same as clusters_track, but each point also has a radius from the CHT.
    Clusters can only be merged if their radii (plus max_dist slack) overlap."""
    # [cx, cy, radius, set_of_frames]  center and radius follow the most recent point
    clusters = []
    for x, y, frame, *_, radius in sorted(points, key=lambda p: p.frame):
        best = None
        best_dist = None
        for cluster in clusters:
            cx, cy, cradius, _ = cluster
            reach = cradius + radius + max_dist
            dist = (x - cx) ** 2 + (y - cy) ** 2
            # closest cluster within a threshold that widens with both radii
            if dist <= reach**2 and (best is None or dist < best_dist):
                best = cluster
                best_dist = dist

        if best is not None:
            # extend the track and move its center and radius to the latest detection
            best[0], best[1], best[2] = x, y, radius
            best[3].add(frame)
        else:
            clusters.append([x, y, radius, {frame}])

    return sum(1 for c in clusters if len(c[3]) >= min_frames)

# Sklearn density clustering on (x, y) detections, with DBSCAN or HDBSCAN; then
# count clusters that persist for a minimum number of frames.
def _count_persistent(points, labels, min_frames):
    """Count clusters (label -1 is noise) that span at least min_frames distinct frames."""
    bubbles = 0
    for label in set(labels):
        if label == -1:
            continue
        frames = {points[i].frame for i, l in enumerate(labels) if l == label}
        if len(frames) >= min_frames:
            bubbles += 1
    return bubbles

def clusters_dbscan(points, max_dist=15, min_frames=3):
    """Density clustering: a bubble is a dense blob of (x, y) detections, stray detections
    become noise. eps = max_dist, then keep blobs that persist min_frames frames."""
    if len(points) < 2:  # a lone detection is noise
        return 0
    xy = np.array([(p.x, p.y) for p in points], dtype=float)
    labels = DBSCAN(eps=max_dist, min_samples=2).fit_predict(xy)
    return _count_persistent(points, labels, min_frames)

def clusters_hdbscan(points, max_dist=15, min_frames=3):
    """Like DBSCAN but picks the density scale itself(no eps guess). max_dist
    only merges clusters closer than that. Keeps blobs that persist min_frames frames."""
    if len(points) < 2:  # HDBSCAN needs >=2 samples; a lone detection is noise
        return 0
    xy = np.array([(p.x, p.y) for p in points], dtype=float)
    try:
        labels = HDBSCAN(min_cluster_size=2, cluster_selection_epsilon=float(max_dist),
                         copy=True).fit_predict(xy)
    except Exception:
        # sklearn HDBSCAN can crash on too few points. treat that as no resolvable bubble rather
        # than killing whole validation run.
        return 0
    return _count_persistent(points, labels, min_frames)

def nbub(ev, algo=clusters_greedy):
    """One event's detections (dict of NumPy columns) -> estimated bubble count."""
    if len(ev.get("cam", [])) == 0:
        return 0

    cam = ev["cam"]
    pos = ev["pos"]
    frame = ev["frame"]
    significance = ev["significance"]
    radius = ev["radius"]
    counts = []

    # count bubbles per camera, trust the cam that saw the most
    for c in (1, 2, 3):
        mask = cam == c
        pts = [Point(x, y, f, s, r)
               for (x, y), f, s, r in zip(pos[mask], frame[mask], significance[mask], radius[mask])]
        counts.append(algo(pts))

    return max(counts) if counts else 0
