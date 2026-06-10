from sbcbinaryformat import Streamer, Writer
from collections import Counter, defaultdict
import numpy as np






def getProjMat(camNum):
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



    '''
    Projection matrices for the three cameras, generated using OpenCV calibration, with
    SiPMs and fiducial markings as training points
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


    if valid_cam_count < 2:
        return np.array([-999,-999, -999])

    A = np.array(A)

    _, _, Vt = np.linalg.svd(A)
    X = Vt[-1]
    X = X / X[3]

    return X[:3]




'''
Check how many bubbles the bubble finder detected in a given event, returns estimated count
'''
def bubble_mult(bubble_data):
    '''
    bubble_data: bubble finder output dictionary

    Returns:
        int of estimated bubble multiplicity
    '''
    # grab info

    def to_frame(x):
        return int(x[0]) if isinstance(x, (list, tuple)) and len(x) > 0 else int(x)

    frames = [to_frame(f) for f in bubble_data["frame"]]
    cams   = [tuple(c) if isinstance(c, list) else c for c in bubble_data["cam"]]

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
    n = 5
    seq = [(f, c) for f, c  in zip(frames, cams) if ( f >= firstFrame and f <= firstFrame + (10 + n)) ]
    if not seq:
        return -1

    mult = Counter(seq)  # {(frame, cam): multiplicity}
    byCamDict = defaultdict(dict)
    lastSeen = set()
    for f, c in seq:
        if (f,c) not in lastSeen:
            byCamDict[c][f] = mult[(f,c)]
            lastSeen.add((f,c)) 
   
    
    sortedByMult = sorted(mult.keys(), key = lambda k: mult[k], reverse=True)
    checked  = []
    for f0, c0 in sortedByMult:
        if ( (f0,c0) in checked):
            continue
        checked.append((f0,c0))
        m0 = mult[(f0,c0)]
        ok = True
        for offset in range(n):    
            if  mult[f0 + offset, c0] < mult[f0, c0]:
                ok = False
                break
        if ok:
            return m0
    return 0



'''
Pull bubble 2D position data from bubble finder. Uses positions from first frame where at
least two cameras are defined. Undefined camera stored as np.nan. Does NOT work
correctly for multi-bubble events yet
'''

def pull_bubble_coords(bubble_data):
    '''
    bubble_data: bubble finder output dictionary

    Returns:
        2D coordinates of bubble, [cam1x cam1y cam2x cam2y cam3x cam3y]
    '''
    run_bubbles = bubble_data

    cams = np.array([c[0] for c in run_bubbles['cam']])
    frames = np.array([f[0] for f in run_bubbles['frame']])
    pos = np.array(run_bubbles['pos'])

    if len(frames) == 0:
        returnList = []
        for i in range(50):
            returnList.append(np.full(6, np.nan)) 
            # No found bubbles
            return returnList

    # In order of frames
    frames_ordered = np.argsort(frames)
    cams = cams[frames_ordered]
    pos = pos[frames_ordered]
    frames = frames[frames_ordered]
    
    unique_frames = np.unique(frames)
    coordsToReturn  = []
    for frame in range(50):
        if frame in unique_frames:

            pick_frame = (frames == frame)

            cams_f = cams[pick_frame]
            pos_f  = pos[pick_frame]

            # Need at least 2 cams
            if len(np.unique(cams_f)) < 2:
                output = np.full(6, np.nan)
                for i in range(6):
                    output[i] = -999
                coordsToReturn.append(output)
                continue
        
            output = np.full(6, np.nan)
            used_cams = set()

            #Fill available cameras
            for cam_id, (x, y) in zip(cams_f, pos_f):
                if cam_id in used_cams:
                # multi-bubble?
                    continue

                used_cams.add(cam_id)

                if cam_id == 1:
                    output[0:2] = [x, y]
                elif cam_id == 2:
                    output[2:4] = [x, y]
                elif cam_id == 3:
                    output[4:6] = [x, y]

            coordsToReturn.append(output)
        else:
            coordsToReturn.append(np.full(6,np.nan))
    return coordsToReturn



def reconstruct_2D_to_3D(data):

    """
    OUTPUT KEY:
    -999, -999, -999: Only a single cam is defined
    -1000,-1000,-1000: Possible multi bubble event
    NaN, NaN, NaN: No bubble finder info for this frame/event
    """

    """
    bubble mult check:
        if failed, return list of len 50 of lists of 3 -999 values
    """

    if "bubble" in data["analysis"]:
        bubble_data = data["analysis"]["bubble"]
        
        if bubble_mult(bubble_data) != 1:
            coordsToReturn = []
            for i in range(50):
                coordsToReturn.append(np.full(3, -999))
            return {"coords_3D": coordsToReturn}
        
        coordsToReturn = []
        """
        if frame not in bubble finder
            add NaN
            continue
        grab 2d coords
        get 3d coord
        add to list
        """
        # Pull 2D coordinate
        coords_2D = pull_bubble_coords(bubble_data)
        for coord in coords_2D:
            # Reconstruct
            if coord[0] == -999 or np.isnan(coord).any():
                coordsToReturn.append(np.full(3,np.nan))
                continue
            coords_3D = triangulate_multi_cam_LS(coord)
            print(coord)
            print(coords_3D)
            coordsToReturn.append(coords_3D)

        return {"coords_3D": coordsToReturn}

    else:
        coordsToReturn = []
        for i in range(50):
            coordsToReturn.append(np.full(3,np.nan))

        print(len(coordsToReturn))
        return {"coords_3D": coordsToReturn}



