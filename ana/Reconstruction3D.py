from sbcbinaryformat import Streamer, Writer
from collections import Counter, defaultdict
import numpy as np
import os

'''
Projection matricies for each camera, generated using OpenCV calibration, with SiPMs and fiducial markings as training points

'''
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
Check how many bubbles the bubble finder detected in a given event, returns estimated count. 0 if no bubbles found
'''
def bubble_mult(bubble_data, frameCount):
    '''
    bubble_data: bubble finder output dictionary

    Returns:
        int of estimated bubble multiplicity
    '''
    # grab info
    
    frames = [int(f[0]) if isinstance(f, list) else int(f) for f in bubble_data["frame"]]
    cams   = [int(c[0]) if isinstance(c, list) else int(c) for c in bubble_data["cam"]]
    sigs   = [float(s[0])if isinstance(s,list) else float(s) for s in bubble_data["significance"]]
    
    # find first mutli cam frame
    firstFrame = -1
    idx = sorted(range(len(frames)), key=lambda i: frames[i])
    seen = set()
    for i in idx:
        seen.add(cams[i])
        if len(seen) >= 2:
            firstFrame = frames[i]
            break

    # get all camera frame pairs within a range of the first mutli cam event
    n = 2
    seq = [(f, c) for f, c, s  in zip(frames, cams, sigs) if ( f >= firstFrame and f <= firstFrame + (10 + n) and s >= 0.75) ]
    if not seq:
        return -1

    mult = Counter(seq)  # {(frame, cam): multiplicity}
    
    sortedByMult = sorted(mult.keys(), key = lambda k: (-mult[k],k))
    checked  = []
    m0 = 0
    for f0, c0 in sortedByMult:
        if ( (f0,c0) in checked):
            continue
        checked.append((f0,c0))
        m0 = mult[(f0,c0)]
        ok = True
        for offset in range(n):    
            if f0 + offset >= frameCount:
                break
            disagree = mult[f0 + offset, 1] < m0
            disagree+= mult[f0 + offset, 2] < m0 
            disagree+= mult[f0 + offset, 3] < m0
            if disagree >= 2:
                ok = False
                break
        if ok:
            return m0
    return m0


'''
Pull bubble 2D position data from bubble finder. Uses positions from first frame where at
least two cameras are defined. Undefined camera stored as np.nan. Does NOT work
correctly for multi-bubble events yet
'''

def pull_bubble_coords(bubble_data, frameCount):
    '''
    bubble_data: bubble finder output dictionary

    Returns:
        List of 2D coordinates of bubble for each frame, [cam1x cam1y cam2x cam2y cam3x cam3y]
    '''
    run_bubbles = bubble_data

    cams = np.array([int(c[0]) for c in run_bubbles['cam']])
    frames = np.array([int(f[0]) for f in run_bubbles['frame']])
    sigs = np.array([float(s[0]) for s in run_bubbles['significance']])
    pos = np.array(run_bubbles['pos'])

    if len(frames) == 0:
        returnList = []
        for i in range(frameCount):
            returnList.append(np.full(6, np.nan), i) 
            # No found bubbles
        return returnList

    # In order of frames
    frames_ordered = np.argsort(frames)
    cams = cams[frames_ordered]
    pos = pos[frames_ordered]
    sigs = sigs[frames_ordered]
    frames = frames[frames_ordered]
    
    unique_frames = np.unique(frames)
    coordsToReturn  = []
    maxFrame = frameCount
    for frame in range(0,maxFrame):
        if frame in unique_frames:

            pick_frame = (frames == frame)

            cams_f = cams[pick_frame]
            pos_f  = pos[pick_frame]
            sig = sigs[pick_frame]
            # Need at least 2 cams
            if len(np.unique(cams_f)) < 2:
                output = np.full(6, np.nan)
                coordsToReturn.append((output,frame))
                continue
            

            output = np.full(6, np.nan)
            used_cams = set()
            used_sigs = []
            
            #Fill available cameras
            for cam_id, (x, y), s in zip(cams_f, pos_f, sigs):
                if cam_id in used_cams:
                    #multi-bubble?
                    if not s > used_sigs[-1]:
                        continue
            
                used_cams.add(cam_id)
                used_sigs.append(s)
                if cam_id == 1:
                    output[0:2] = [x, y]
                elif cam_id == 2:
                    output[2:4] = [x, y]
                elif cam_id == 3:
                    output[4:6] = [x, y]
            coordsToReturn.append((output, frame))
        else:
            coordsToReturn.append((np.full(6,np.nan),frame))
    return coordsToReturn


def reproj(P,x):
    x = x/25.4
    X_h = np.append(x,1.0)
    proj = P @ X_h
    proj = proj[:2]/ proj[2]
    return proj



def reconstruct_2D_to_3D(data):
    frameCount = 0
    for maybekey in (data["cam"]["c1"].keys()):
        if "frame" in maybekey:
            frameCount += 1
    # checking if the bubble finder ran
    if "bubble" in data["analysis"]:
        bubble_data = data["analysis"]["bubble"]
        # if it ran and the estimated multiplicity isnt 1, this could be a multi bubble so we ignore it.
        mult = bubble_mult(bubble_data, frameCount)
        if mult > 1:
            coordsToReturn = []
            coordsToReturn.append((-1000,-1000,-1000))
            frames = []
            reprojErrors = []
            for i in range(0,frameCount):
                coordsToReturn.append((-1000,-1000,-1000))
                reprojErrors.append(np.nan)
                frames.append(i)
            return {"coords_3D": coordsToReturn, "frame": frames, "reprojError": reprojErrors}
        # list of 3d coordinates to return to event dealer
        coordsToReturn = []
    
        # Pulls all 2D coordinates
        coords_2D  = pull_bubble_coords(bubble_data, frameCount)
        frames = []
        reprojErrors = []
        # for every frame there is a set of 2d coordinates, each one corresponding to a certian cameras bubble location
        for coord in coords_2D:
            # if the camera didnt have a bubble, we should just ignore this frame and more on
            if isinstance(coord, int) or len(coord) != 2:
                coordsToReturn.append(np.full(3,np.nan))
                frames.append(len(frames))
                reprojErrors.append(np.nan)
                continue
            frames.append(coord[1])
            nancheck = 0
            for i in coord[0]:
                if np.isnan(i):
                    nancheck += 1
            if len(coord) != 2 or len(coord[0]) != 6 or  coord[0][0] <= -999 or nancheck >= 4:
                coordsToReturn.append((-999,-999,-999))
                reprojErrors.append(np.nan)
                continue
            # triangulate the bubble into 3d space, then add it to the list to return
            coords_3D = triangulate_multi_cam_LS(coord[0])
            coordsToReturn.append(coords_3D)
            reprojCoord = reproj(getProjMat(1), coords_3D)
            reprojError = 0
            count = 0
            if not np.isnan(coord[0][0]) and not np.isnan(coord[0][1]):
                reprojError += np.linalg.norm(reproj(getProjMat(1), coords_3D) - (coord[0][0], coord[0][1]))
                count += 1
            if not np.isnan(coord[0][2]) and not np.isnan(coord[0][3]):
                reprojError += np.linalg.norm(reproj(getProjMat(2), coords_3D) - (coord[0][2], coord[0][3]))
            if not np.isnan(coord[0][4]) and not np.isnan(coord[0][5]):
                reprojError += np.linalg.norm(reproj(getProjMat(3), coords_3D) - (coord[0][4], coord[0][5]))
            if count != 0:
                reprojError /= count
            else:
                reprojError = np.nan
            reprojErrors.append(1/reprojError)
        

        return {"coords_3D": coordsToReturn, "frame": frames, "reprojError": reprojErrors}

    else:
        coordsToReturn = []
        coordsToReturn.append(np.full(3,np.nan))
        frames = []
        reprojErrors = []
        for i in range(0,frameCount):
            coordsToReturn.append(np.full(3,np.nan))
            frames.append(i)
            reprojErrors.append(np.nan)
        return {"coords_3D": coordsToReturn, "frame": frames, "reprojError": reprojErrors}

