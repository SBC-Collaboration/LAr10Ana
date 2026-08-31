from sbcbinaryformat import Streamer, Writer
import numpy as np
import matplotlib.pyplot as plt
from skimage.draw import circle_perimeter, disk
from skimage.measure import label, regionprops
import diplib as dip
from skimage.restoration import rolling_ball

"""
Args:
  ev: event
  cam: camera
  noise_thresh: diff values below this threshold will be set to zero

Returns:
  bub_dict: dictionary of lists, where each row is a bubble frame from one camera
    cam (int): camera number of this bubble
    pos (float, 2): x and y axis of the pixel position of the bubble
    radius (float): radius of the bubble in pixels
    frame (int): frame number of this bubble
    second_good_cand (int)(but only 0 or 1): 1 if the second most prominent circular feature falls 
    near a sufficiently bright and large region of connected pixels; indicates likely at 
    least one more bubble in frame; can be 100% ignored by analyzers, only used to get confidence
    confidence (float): for each camera, the percentage of frames where a bubble has been found and
    there is not likely a second bubble or prominent circular structure in the event; for handscanned,
    confirmed single bubble events, low confidence indicates a particularly noisy image and/or the 
    likelihood that the real bubble has been missed.  Cut for low confidence should be about ~70% 
"""


out_keys = ['rad','pos','frame','cam','second_good_cand']
   
def _new_bub_dict():
    return dict([(key, []) for key in out_keys])
    
def FindBubbles(ev, cam, noise_thresh, bub_dict=None):

    if not ev['cam'][f'c{cam}']['loaded']:
        return bub_dict

    stable_pos = 0
    prev_pos = [0,0]
    constrained_search = False
    prev_cent = [0,0]
    bubless_frames = 0
    bad = False
   
    #get mask for bubble region based on camera
    refIm = np.float32(np.average(ev['cam'][f'c{cam}']['frame0'],axis=2))
    imShape = refIm.shape
    Tshape = imShape[::-1]
        
    if cam==1:
        circy, circx = disk((375, 595), 290, shape=imShape)
        coord_mask = 1.9*circx-1200<circy
        mask_circx = circx[coord_mask]
        mask_circy = circy[coord_mask]
    elif cam==2:
        mask_circy, mask_circx = disk((420, 670), 290, shape=imShape)
    elif cam==3:
        mask_circy, mask_circx = disk((440, 690), 290, shape=imShape)

    if bub_dict is None:
        bub_dict = _new_bub_dict()
    elif not isinstance(bub_dict, dict):
        raise ValueError("bub_dict must be a dictionary")
    elif len(bub_dict) == 0:
        bub_dict = _new_bub_dict()
    elif not all(key in bub_dict for key in out_keys):
        raise ValueError("bub_dict does not contain all required keys: %s" % out_keys)
    elif not all(isinstance(bub_dict[key], list) for key in out_keys):
        raise ValueError("All values in bub_dict must be lists")
    elif not all(len(bub_dict[key]) == len(bub_dict["rad"]) for key in out_keys):
        raise ValueError("All lists in bub_dict must have the same length as bub_dict['rad']")

    keys = ev['cam'][f'c{cam}'].keys()
    frames = 0
    for key in keys:
        if 'frame' in key:
            frames+=1

    for i in range(1,frames):
        
        im_num = frames-i

        diff = np.zeros((imShape[0],imShape[1]))
        if constrained_search==False:
            
            prevIm = np.float32(np.average(ev['cam'][f'c{cam}'][f'frame{im_num-1}'],axis=2))
            thisIm = np.float32(np.average(ev['cam'][f'c{cam}'][f'frame{im_num}'],axis=2))
            
            preMask_prevDiff = abs(thisIm-prevIm)
            preMask_refDiff = abs(thisIm-refIm)
            
            preMask_prevDiff[preMask_prevDiff<noise_thresh] = 0
            preMask_refDiff[preMask_refDiff<noise_thresh] = 0
            
            preMask_prevDiff-=dip.GetSinglePixels(preMask_prevDiff > 0)
            preMask_refDiff-=dip.GetSinglePixels(preMask_refDiff > 0)

            if np.std(preMask_prevDiff)<np.std(preMask_refDiff):
                diff[mask_circy,mask_circx] = preMask_prevDiff[mask_circy,mask_circx]  
            else:
                diff[mask_circy,mask_circx] = preMask_refDiff[mask_circy,mask_circx]
                
        else:
            
            thisIm = np.float32(np.average(ev['cam'][f'c{cam}'][f'frame{im_num}'],axis=2))
            preMask_refDiff = abs(thisIm-refIm)
            preMask_refDiff[preMask_refDiff<noise_thresh] = 0
            preMask_refDiff-=dip.GetSinglePixels(preMask_refDiff > 0)
            
            diff[con_y,con_x] = preMask_refDiff[con_y,con_x]
            
        diff-=dip.GetSinglePixels(diff > 0) 
        
        filt = rolling_ball(diff, radius = 1)
        filt[filt<np.mean(filt)+3*np.std(filt)] = 0
        labelIm = filt
    
        #connectivity = 2 allows pixels to count as connected if they are diagonal from each other
        labeled = label(labelIm>0, connectivity = 2)
        
        #get properties of labeled regions
        props = regionprops(labeled, intensity_image = labelIm)
        if len(props) == 0:
            continue
            
        #get region properties: area, mean pixel intensity, length
        intensities = np.array([prop.intensity_mean for prop in props])
        areas = np.array([prop.area for prop in props])
        lengths = np.array([prop.axis_major_length for prop in props])

        if constrained_search==False:
            
            intensity_thresh = np.mean(intensities)+3*np.std(intensities)
            area_thresh = np.mean(areas)+3*np.std(areas)
            length_thresh = np.mean(lengths)+3*np.std(lengths)

            intensities[intensities<intensity_thresh] = 0
            areas[areas<area_thresh] = 0
            lengths[lengths<length_thresh] = 0
            
            #don't want to stabilize position and search around poor candidates
            if len(intensities[intensities>0])==0 or len(areas[areas>0])==0 or len(lengths[lengths>0])==0:
                continue
                
        elif np.max(areas)<10: #end backwards scan once region areas get this small, bubs usually appear at 6-10ish pixels
            return bub_dict

        #as long as constrained search is true, always assume a bubble somewhere in the search region unless the area
        #is too small or we've reached 3 consecutive frames below threshold

        #choose highest scoring candidate
        scores = (areas/np.max(areas) + intensities/np.max(intensities) + lengths/np.max(lengths))/3
        cand_idx = np.argmax(scores)
        
        cand_region = props[cand_idx]
        score_thresh = np.mean(scores) + 3*np.std(scores)
        cand_regions = np.array(props)[scores>=score_thresh]

        #estimate rad cands
        min_est_rad = np.round(cand_region.axis_major_length/2)
        if min_est_rad - 2 <= 2:
            min_rad = 2
            max_rad = 6
        else:
            min_rad = min_est_rad - 2
            max_rad = min_est_rad + 3
        rad_cands = np.arange(min_rad, max_rad,1)
        
        #perform CHT
        for rad in rad_cands:
            circx, circy = circle_perimeter(600, 400, int(rad), shape=Tshape)
            dx = circx-600
            dy = circy-400
            offsets = [(dx[i],dy[i]) for i in range(len(dx))]
        
            this_layer = np.zeros((imShape[0], imShape[1]))
            for offset in offsets:
                this_layer += np.roll(labelIm, offset,(0,1))
        
            if rad==rad_cands[0]:
                accum = this_layer
            else:
                accum = np.dstack((accum, this_layer))
        
        accum_shape = accum.shape
        
        rcy, rcx = cand_region.centroid
        ry, rx = disk((rcy,rcx), 20, shape=imShape)
        regIm = np.zeros(accum_shape)
        regIm[ry,rx] = accum[ry, rx]
    
        pcy, pcx, rad_ind = np.unravel_index(np.argmax(regIm), accum_shape)
        prad = rad_ind + rad_cands[0]

        buby, bubx = disk((pcy,pcx), prad, shape=imShape)

        if constrained_search==True and np.mean(filt[buby,bubx])<np.mean(filt)+2.5*np.std(filt): 
            bubless_frames+=1
        elif constrained_search==True and np.mean(filt[buby,bubx])>=np.mean(filt)+2.5*np.std(filt): 
            bubless_frames = 0 #stop search after 3 consecutive frames below threshold 

        #delete this bubble from the accumulator array and see how the next highest peak compares
        peak_votes = np.max(regIm)
        bub_neighborhood_y, bub_neighborhood_x = disk((pcy,pcx), prad+20, shape=imShape)
        accum[bub_neighborhood_y, bub_neighborhood_x] = 0
        next_peak_votes = np.max(accum)
        vote_rat = next_peak_votes/peak_votes

        #see if next accum peak falls within any of the other high-scoring candidate regions
        npcy, npcx, _ = np.unravel_index(np.argmax(accum), accum_shape)
        second_good_cand = 0
        for reg in cand_regions:
            cand_cent = reg.centroid
            dist = np.sqrt((npcx-cand_cent[1])**2+(npcy-cand_cent[0])**2)
            if dist<20:
                second_good_cand = 1
                break

        if bubless_frames==3:
            return bub_dict
        elif constrained_search==True and vote_rat<1:
            bub_dict["rad"].append([prad])
            bub_dict["cam"].append([cam])
            bub_dict["frame"].append([im_num])
            bub_dict["pos"].append([pcx,pcy])
            bub_dict['second_good_cand'].append([second_good_cand])

        if constrained_search==False:
            cent = cand_region.centroid
            if np.sqrt((cent[0]-prev_cent[0])**2 + (cent[1]-prev_cent[1])**2)<20:
                if vote_rat<1:
                    stable_pos+=1
                else:
                    stable_pos = 0
            if vote_rat<1:
                prev_pos = [pcx,pcy]
                prev_cent = cent
            if stable_pos==2: #consistent between images twice --> 3 consecutive consistent images
                constrained_search = True
                con_x, con_y = disk((pcx, pcy), prad+10, shape=Tshape)
        else:
            con_x, con_y = disk((pcx, pcy), prad+10, shape=Tshape) 

    return bub_dict

def BubbleFinder(ev, noise_thresh = 5):
    
    out = _new_bub_dict()
    out = FindBubbles(ev, 1, noise_thresh, bub_dict=out)
    out = FindBubbles(ev, 2, noise_thresh, bub_dict=out)
    out = FindBubbles(ev, 3, noise_thresh, bub_dict=out)

    if len(out["rad"]) == 0:
        raise ValueError("No bubbles found in event")

    out['confidence'] = []

    #calculate and add confidence score for each camera;
    #confidence that there are no more bubbles visible in each set of camera frames
    cam = np.array(out['cam']).ravel()
    for c in range(1,4):
        
        cMask = cam==c
        if len(cam[cMask])==0: #no bubs found for this camera
            continue
    
        #1 - how often the next most prominent circular feature lines up with a large 
        #and bright region of connected pixels, or how often there is likely at least 
        #one more bubble to be found in the image
        sgc = np.array(out['second_good_cand']).ravel()[cMask]
        conf = 1 - len(sgc[sgc>0])/len(sgc) 
        
        for i in range(len(cam[cMask])):
            out['confidence'].append([conf])

    #use confidence scores to get t0 frame and cams containing t0 frame
    out['t0_frame'] = []
    out['t0_cams'] = []
    
    c1Mask = cam==1
    c2Mask = cam==2
    c3Mask = cam==3
    confs = np.array(out['confidence']).ravel()
    rs = np.array(out['rad']).ravel()
    fs = np.array(out['frame']).ravel()
    
    conf1 = confs[c1Mask]
    conf2 = confs[c2Mask]
    conf3 = confs[c3Mask]
    
    #only want to consider cams with confidence score above 70%
    if len(conf1)>0 and conf1[0]>0.7:
        r1 = rs[c1Mask]
        f1 = fs[c1Mask]
    else:
        r1 = [0]
        f1 = [100]
    
    if len(conf2)>0 and conf2[0]>0.7:
        r2 = rs[c2Mask]
        f2 = fs[c2Mask]
    else:
        r2 = [0]
        f2 = [100]
    
    if len(conf3)>0 and conf3[0]>0.7:
        r3 = rs[c3Mask]
        f3 = fs[c3Mask]
    else:
        r3 = [0]
        f3 = [100]
    
    #get earliest frame where one camera dips to the lowest radius it finds
    earliest = np.min([f1[np.argmin(r1)],f2[np.argmin(r2)],f3[np.argmin(r3)]])
    #get closest frame to this where at least 2 cams have a detection
    f12 = np.intersect1d(f1,f2)
    f23 = np.intersect1d(f2,f3) 
    f13 = np.intersect1d(f1,f3)
    
    if len(f12)==0:
        f12 = [100]
    if len(f23)==0:
        f23 = [100]
    if len(f13)==0:
        f13 = [100]
        
    f12_earliest = f12[np.argmin(abs(f12 - earliest))]
    f23_earliest = f23[np.argmin(abs(f23 - earliest))]
    f13_earliest = f13[np.argmin(abs(f13 - earliest))]
    
    t0_frame = np.min([f12_earliest,f23_earliest,f13_earliest])
    ct0 = np.argmin([f12_earliest,f23_earliest,f13_earliest])
    if ct0==0:
        t0_cams = [1,2]
    elif ct0==1:
        t0_cams = [2,3]
    elif ct0==2:
        t0_cams = [1,3]

    for i in range(len(cam)):
        out['t0_frame'].append([t0_frame])
        out['t0_cams'].append(t0_cams)

    return out
