from sbcbinaryformat import Streamer, Writer
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



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

    P1 = np.array([[-1.05302109e+02, -7.02444185e+02, -3.34577970e+02,  5.72535995e+03],
     [-5.51213766e+02,  2.58404210e+01, -3.45420423e+02,  3.46877200e+03],
     [ 5.46200003e-02, -3.31725499e-01, -9.41793422e-01,  8.93247437e+00]])
    P2 = np.array([[ 6.24551374e+02,  2.05426176e+02, -4.20327029e+02,  6.08648299e+03],
     [ 2.38142395e+02, -5.16479247e+02, -3.98154885e+02,  3.57897785e+03],
     [ 1.75014306e-01,  8.21425255e-02, -9.81133323e-01,  8.45879059e+00]])
    P3 = np.array([[-4.46470566e+02,  4.77173422e+02, -4.42541834e+02,  5.80637791e+03],
     [ 3.67166284e+02,  4.75216795e+02, -4.43358757e+02,  3.38952193e+03],
     [-9.35610736e-02,  1.48157021e-01, -9.84528223e-01,  7.62754209e+00]])

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
        return np.array([np.nan, np.nan, np.nan])

    A = np.array(A)

    _, _, Vt = np.linalg.svd(A)
    X = Vt[-1]
    X = X / X[3]

    return X[:3]






'''
Pull bubble 2D position data from bubble finder. Uses positions from first frame where at
least two cameras are defined. Undefined camera stored as np.nan. Does NOT work
correctly for multi-bubble events yet
'''

def pull_bubble_coords(bubble_data,event_number):
    '''
    bubble_data: bubble finder output dictionary
    event_number: event number in the run

    Returns:
        2D coordinates of bubble, [cam1x cam1y cam2x cam2y cam3x cam3y]
    '''
    run_bubbles = bubble_data
    run_bubbles = run_bubbles.to_dict()

    pick_event = run_bubbles['ev'] == event_number

    if not np.any(event_number): # Event number not found
        return np.full(6, np.nan) 

    cams = run_bubbles['cam'][pick_event]
    pos = run_bubbles['pos'][pick_event]
    frames = run_bubbles['frame'][pick_event]


    # In order of frames
    frames_ordered = np.argsort(frames)
    cams = cams[frames_ordered]
    pos = pos[frames_ordered]
    frames = frames[frames_ordered]
    
    unique_frames = np.unique(frames)

    for frame in unique_frames:

        pick_frame = frames == frame

        cams_f = cams[pick_frame]
        pos_f  = pos[pick_frame]

        # Need at least 2 cams
        if len(np.unique(cams_f)) < 2:
            continue
        
        output = np.full(6, np.nan)

        # Fill available cameras
        for cam_id, (x, y) in zip(cams_f, pos_f):

            if cam_id == 1:
                output[0] = x
                output[1] = y

            elif cam_id == 2:
                output[2] = x
                output[3] = y

            elif cam_id == 3:
                output[4] = x
                output[5] = y

        return output

    # If no frame with 2+ cams
    print('3D reconstruction failed: no frame with 2+ cams')
    return np.full(6, np.nan)





def reconstruct_2D_to_3D(data):
    
    if "bubble" in data["analysis"].keys():
        bubble_data = data["analysis"]["bubble"]
        event = data['event_info']['event_id'][0]
        
        # Pull 2D coordinate
        coords_2D = pull_bubble_coords(bubble_data,event_number)
        # Reconstruct
        coords_3D = triangulate_multi_cam_LS(coords_2D)

        return coords_3D

    else:
        print('3D Reconstruction failed: bubble finder data not found')
        return np.full(3, np.nan)



