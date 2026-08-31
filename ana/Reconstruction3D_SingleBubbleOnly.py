from sbcbinaryformat import Streamer, Writer
import numpy as np
import os

'''
3D reconstruction for the single-bubble-only bubble finder.

Expected bubble finder output dictionary (data["analysis"]["bubble"]):
    rad               (N,1)  bubble radius in pixels
    pos               (N,2)  [x, y] pixel position of the bubble
    frame             (N,1)  frame number of this row
    cam               (N,1)  camera number of this row
    second_good_cand  (N,1)  internal to the finder, ignored here
    confidence        (N,1)  per-camera confidence, constant across a camera's rows.
                             Do NOT trust a camera at ~0.70 or below.
    t0_frame                 the t0 frame. 100 means no bubbles were found with
                             high confidence anywhere in the event.
    t0_cams                  the two cameras that have the t0 frame with high confidence.

Note the rename from the old finder: 'radius' -> 'rad', and 'significance' is gone,
replaced by the per-camera 'confidence'. The multiplicity estimate is also gone: this
finder is single-bubble-only, so there is nothing to disambiguate.
'''

# a camera at or below this confidence is not trusted
CONFIDENCE_THRESHOLD = 0.70

# t0_frame sentinel meaning "no bubbles found with high confidence"
NO_T0_FRAME = 100

# analysis keys to look for the single bubble finder output under, in order
BUBBLE_KEYS = ("bubble_single", "bubble")

# coords_3D sentinels
NO_HIGH_CONF_T0 = -1000.0   # event has no high-confidence t0 at all
TOO_FEW_CAMS = -999.0       # this frame has bubbles but < 2 trusted cameras


def _distance_to_wall(x, y, z):
    '''
    x,y,z: 3D coordinates in mm

    Returns:
        distance to wall in mm
    '''

    in_to_mm = 25.4
    r_inside = 4.525 * in_to_mm
    thickness = 0.2 * in_to_mm
    r_outside = r_inside + thickness
    z_low = -12*in_to_mm
    z_high = (14.72-15.358)*in_to_mm

    angle_neck = 1.19367  # angle where the neck and dome meet
    r_neck = 2*in_to_mm
    r_neck_inside = r_neck - thickness
    x_neck, z_neck = 2.725*in_to_mm, z_high  # center of neck circle
    r_dome = 9.4*in_to_mm
    r_dome_inside = r_dome - thickness
    z_dome = (7.84 - 15.358)*in_to_mm  # center of dome circle
    R = r_inside  # for r^2/R

    if np.isnan(x) or np.isnan(y) or np.isnan(z):
        return -1

    r = np.sqrt(x**2 + y**2)
    angle = np.arctan2(z - z_dome, r)  # angle based on dome center
    if z <= z_high:
        return r_inside - r
    elif 0 <= angle and angle <= angle_neck:
        return r_neck_inside - np.sqrt((r - x_neck)**2 + (z - z_neck)**2)
    elif angle > angle_neck and angle <= np.pi/2:
        return r_dome_inside - np.sqrt(r**2 + (z - z_dome)**2)
    else:
        return -1

def _distance_to_wall_arr(arr):
    '''
    arr: 3D coordinates in mm, shape (N,3)

    Returns:
        distance to wall in mm, shape (N,)
    '''
    return np.array([_distance_to_wall(*pos) for pos in arr])

def getProjMat(camNum):
    '''
    camNum: 1,2, or3
    Returns:
        4x3 matrix or np.nan if invalid camera
    '''
    if camNum == 1:
        return np.array([[-1.05302109e+02, -7.02444185e+02, -3.34577970e+02,  5.72535995e+03],
                [-5.51213766e+02,  2.58404210e+01, -3.45420423e+02,  3.46877200e+03],
                [ 5.46200003e-02, -3.31725499e-01, -9.41793422e-01,  8.93247437e+00]])


    if camNum == 2:
        return np.array([[ 6.24551374e+02,  2.05426176e+02, -4.20327029e+02,  6.08648299e+03],
                [ 2.38142395e+02, -5.16479247e+02, -3.98154885e+02,  3.57897785e+03],
                [ 1.75014306e-01,  8.21425255e-02, -9.81133323e-01,  8.45879059e+00]])


    if camNum == 3:
        return np.array([[-4.46470566e+02,  4.77173422e+02, -4.42541834e+02,  5.80637791e+03],
                [ 3.67166284e+02,  4.75216795e+02, -4.43358757e+02,  3.38952193e+03],
                [-9.35610736e-02,  1.48157021e-01, -9.84528223e-01,  7.62754209e+00]])
    return np.array([[np.nan,  np.nan, np.nan,  np.nan],
     [np.nan, np.nan, np.nan, np.nan],
     [np.nan,np.nan,np.nan,np.nan]])


'''
Least squares triangulation, 2D to 3D points. Needs 2+ cams
'''

def triangulate_multi_cam_LS(pixel_coords):
    '''
    pixel_coords: [cam1x,cam1y,cam2x,cam2y,cam3x,cam3y] with np.nan where missing cam

    Returns:
        3D point (X,Y,Z) or np.nan if not enough defined points
    '''

    P1 = getProjMat(1)
    P2 = getProjMat(2)
    P3 = getProjMat(3)

    P_mats = [P1, P2, P3]

    pixel_coords = np.asarray(pixel_coords).reshape(3, 2)
    A = []
    valid_cam_count = 0

    for P, (x, y) in zip(P_mats, pixel_coords):

        # Skip camera if either coordinate is np.nan
        if np.isnan(x) or np.isnan(y):
            continue
        valid_cam_count += 1

        A.append(x * P[2] - P[0])
        A.append(y * P[2] - P[1])

    # if there isnt 2 or more cameras, we cant triangulate
    if valid_cam_count < 2:
        return np.array([np.nan,np.nan, np.nan])

    A = np.array(A)

    _, _, Vt = np.linalg.svd(A)
    X = Vt[-1]
    X = X / X[3]
    return X[:3] * 25.4


'''
Helpers for reading the finder output. Every per-row value is written as a
one-element list, but t0_frame / t0_cams may be bare scalars or lists depending
on how the finder was invoked, so unwrap defensively.
'''

def _scalar(value):
    '''
    Unwrap a possibly nested one-element list/array down to a scalar.
    '''
    while isinstance(value, (list, tuple, np.ndarray)):
        if len(value) == 0:
            return None
        value = value[0]
    return value

def _flat_list(value):
    '''
    Flatten a possibly nested list/array into a flat python list.
    '''
    if value is None:
        return []
    if not isinstance(value, (list, tuple, np.ndarray)):
        return [value]
    out = []
    for item in value:
        out.extend(_flat_list(item))
    return out


def get_t0_frame(bubble_data):
    '''
    bubble_data: bubble finder output dictionary

    Returns:
        int t0 frame, or NO_T0_FRAME (100) if no high confidence bubbles
    '''
    t0 = _scalar(bubble_data.get("t0_frame"))
    if t0 is None:
        return NO_T0_FRAME
    return int(t0)


def get_t0_cams(bubble_data):
    '''
    bubble_data: bubble finder output dictionary

    Returns:
        sorted list of the cameras that have the t0 frame with high confidence
    '''
    return sorted({int(c) for c in _flat_list(bubble_data.get("t0_cams"))})


def camera_confidence(bubble_data):
    '''
    bubble_data: bubble finder output dictionary

    Returns:
        dict of {cam number: confidence}. Confidence is constant across a
        camera's rows, so the first row for each camera is enough.
    '''
    cams = [int(c) for c in _flat_list(bubble_data.get("cam"))]
    confs = [float(s) for s in _flat_list(bubble_data.get("confidence"))]

    conf_by_cam = {}
    for cam, conf in zip(cams, confs):
        conf_by_cam.setdefault(cam, conf)
    return conf_by_cam


def trusted_cameras(bubble_data):
    '''
    bubble_data: bubble finder output dictionary

    Returns:
        sorted list of cameras whose output we trust. The t0 cameras are high
        confidence by definition, any other camera has to clear the threshold
        on its own.
    '''
    conf_by_cam = camera_confidence(bubble_data)
    trusted = set(get_t0_cams(bubble_data))
    for cam, conf in conf_by_cam.items():
        if conf > CONFIDENCE_THRESHOLD:
            trusted.add(cam)
    return sorted(trusted)


'''
Pull bubble 2D position data from the single-bubble finder. One bubble per
camera per frame, so there is nothing to disambiguate; rows from untrusted
cameras and frames before t0 are dropped.
'''

def pull_bubble_coords(bubble_data, frameCount, cams_to_use=None, t0_frame=None):
    '''
    bubble_data: bubble finder output dictionary
    frameCount: number of frames in the event
    cams_to_use: iterable of trusted camera numbers, defaults to trusted_cameras()
    t0_frame: frames before this have no bubble, defaults to get_t0_frame()

    Returns:
        List of (coords, frame) where coords is
        [cam1x cam1y cam2x cam2y cam3x cam3y] with np.nan for missing cameras
    '''
    if cams_to_use is None:
        cams_to_use = trusted_cameras(bubble_data)
    if t0_frame is None:
        t0_frame = get_t0_frame(bubble_data)
    cams_to_use = set(int(c) for c in cams_to_use)

    cams = np.array([int(c) for c in _flat_list(bubble_data.get("cam"))])
    frames = np.array([int(f) for f in _flat_list(bubble_data.get("frame"))])
    pos = np.array(bubble_data.get("pos"), dtype=float).reshape(-1, 2)

    if len(frames) == 0:
        return [(np.full(6, np.nan), frame) for frame in range(frameCount)]

    coordsToReturn = []
    for frame in range(frameCount):
        output = np.full(6, np.nan)

        # nothing has happened yet before t0
        if frame < t0_frame:
            coordsToReturn.append((output, frame))
            continue

        pick_frame = (frames == frame)
        cams_f = cams[pick_frame]
        pos_f = pos[pick_frame]

        filled_cams = set()
        for cam_id, (x, y) in zip(cams_f, pos_f):
            if cam_id not in cams_to_use:
                continue
            # single bubble per camera per frame, so ignore any repeat
            if cam_id in filled_cams:
                continue
            filled_cams.add(cam_id)
            if cam_id == 1:
                output[0:2] = [x, y]
            elif cam_id == 2:
                output[2:4] = [x, y]
            elif cam_id == 3:
                output[4:6] = [x, y]

        coordsToReturn.append((output, frame))
    return coordsToReturn


def reproj(P,x):
    x = x/25.4
    X_h = np.append(x,1.0)
    proj = P @ X_h
    proj = proj[:2]/ proj[2]
    return proj


def _flat_output(frameCount, fill):
    '''
    Build an output dictionary where every frame carries the same sentinel.
    '''
    return {"coords_3D": [[fill, fill, fill] for _ in range(frameCount)],
            "frame": [[i] for i in range(frameCount)],
            "reprojError": [[np.nan] for _ in range(frameCount)],
            "d_wall": [[np.nan] for _ in range(frameCount)]}


def reconstruct_2D_to_3D(data):
    def _count_frames(cam_data):
        if not cam_data["loaded"]:
            return 0
        return sum(1 for key in cam_data if key.startswith("frame"))

    frameCount = max([_count_frames(data["cam"][cam]) for cam in ["c1", "c2", "c3"]])

    if frameCount == 0:
        raise ValueError("No frame detected in this event")

    # checking if the single bubble finder ran. It is registered as
    # "bubble_single" in EventDealer; fall back to "bubble" only if that output
    # carries the t0 keys this module needs.
    bubble_data = None
    for key in BUBBLE_KEYS:
        if key in data["analysis"] and "t0_frame" in data["analysis"][key]:
            bubble_data = data["analysis"][key]
            break

    if bubble_data is None:
        return _flat_output(frameCount, np.nan)

    t0_frame = get_t0_frame(bubble_data)

    # t0_frame == 100 means the finder found no bubbles it is confident about
    if t0_frame == NO_T0_FRAME:
        return _flat_output(frameCount, NO_HIGH_CONF_T0)

    cams_to_use = trusted_cameras(bubble_data)

    # the two t0 cameras are high confidence by definition, so this should not
    # happen, but without two trusted cameras there is nothing to triangulate
    if len(cams_to_use) < 2:
        return _flat_output(frameCount, NO_HIGH_CONF_T0)

    # list of 3d coordinates to return to event dealer
    coordsToReturn = []
    frames = []
    reprojErrors = []
    d_wall = []

    # Pulls all 2D coordinates
    coords_2D = pull_bubble_coords(bubble_data, frameCount, cams_to_use, t0_frame)

    # for every frame there is a set of 2d coordinates, each one corresponding to a certian cameras bubble location
    for coord, frame in coords_2D:
        frames.append([frame])

        # need at least 2 of the 6 pixel coordinate pairs defined to triangulate
        defined_cams = sum(1 for i in range(0, 6, 2)
                           if not np.isnan(coord[i]) and not np.isnan(coord[i + 1]))

        # before t0 there is simply no bubble yet
        if frame < t0_frame:
            coordsToReturn.append([np.nan, np.nan, np.nan])
            reprojErrors.append([np.nan])
            d_wall.append([np.nan])
            continue

        if defined_cams < 2:
            coordsToReturn.append([TOO_FEW_CAMS, TOO_FEW_CAMS, TOO_FEW_CAMS])
            reprojErrors.append([np.nan])
            d_wall.append([np.nan])
            continue

        # triangulate the bubble into 3d space, then add it to the list to return
        coords_3D = triangulate_multi_cam_LS(coord)
        coordsToReturn.append([float(coords_3D[0]),
                               float(coords_3D[1]),
                               float(coords_3D[2])])

        reprojError = 0
        count = 0
        for cam_id in (1, 2, 3):
            ix = 2 * (cam_id - 1)
            if np.isnan(coord[ix]) or np.isnan(coord[ix + 1]):
                continue
            reprojError += np.linalg.norm(
                reproj(getProjMat(cam_id), coords_3D) - (coord[ix], coord[ix + 1]))
            count += 1
        if count != 0:
            reprojError /= count
        else:
            reprojError = np.nan
        reprojErrors.append([1/reprojError])

        d_wall.append([_distance_to_wall(*coords_3D)])

    return {"coords_3D": coordsToReturn,
            "frame": frames,
            "reprojError": reprojErrors,
            "d_wall": d_wall}
