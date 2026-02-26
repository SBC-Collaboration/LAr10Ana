from sbcbinaryformat import Streamer, Writer
import numpy as np
import matplotlib.pyplot as plt
from skimage.draw import circle_perimeter, disk
from skimage.measure import label, regionprops
import diplib as dip

"""
Args:
  ev: event
  cam: camera
  num_pix_in_neighborhood: radius around the center of an identified bubble that 
    is zeroed out for all radii in the accumulator array before looking at the next bubble candidate;
    prevents the same bubble from being tagged twice 
  noise_thresh: diff values below this threshold will be set to zero

Returns:
  bub_dict: dictionary of lists, where each row is a bubble frame from one camera
    bub_num (int): index of bubbles. Each bubble in each frame from each camera will have its own index
    cam (int): camera number of this bubble
    pos (float, 2): x and y axis of the pixel position of the bubble
    radius (float): radius of the bubble in pixels
    significance (float):
    frame (int): frame number of this bubble
"""

bub_dict_keys = ["bub_num", "cam", "pos", "radius", "significance", "frame"]

def _new_bub_dict():
    return dict([(key, []) for key in bub_dict_keys])

def FindBubbles(ev, cam, num_pix_in_neighborhood, noise_thresh, bub_dict=None):
    if not ev['cam'][f'c{cam}']['loaded']:
        return bub_dict
    
    #get mask for bubble region based on camera
    refIm = np.float32(np.average(ev['cam'][f'c{cam}']['frame0'],axis=2))
    imShape = refIm.shape
    Tshape = imShape[::-1]
    
    if cam==1:
        circy, circx = disk((375, 590), 280, shape=imShape)
        coord_mask = 2.2*circx-1375<circy
        mask_circx = circx[coord_mask]
        mask_circy = circy[coord_mask]
    elif cam==2:
        mask_circy, mask_circx = disk((425, 670), 270, shape=imShape)
    elif cam==3:
        mask_circy, mask_circx = disk((440, 690), 270, shape=imShape)

    #initialize some variables and the bubble dictionary to be returned at the end
    if bub_dict is None:
        bub_dict = _new_bub_dict()
        start_num = 0
    elif not isinstance(bub_dict, dict):
        raise ValueError("bub_dict must be a dictionary")
    elif len(bub_dict) == 0:
        bub_dict = _new_bub_dict()
        start_num = 0
    elif not all(key in bub_dict for key in bub_dict_keys):
        raise ValueError("bub_dict does not contain all required keys: %s" % bub_dict_keys)
    elif not all(isinstance(bub_dict[key], list) for key in bub_dict_keys):
        raise ValueError("All values in bub_dict must be lists")
    elif not all(len(bub_dict[key]) == len(bub_dict["bub_num"]) for key in bub_dict_keys):
        raise ValueError("All lists in bub_dict must have the same length as bub_dict['bub_num']")
    else:
        start_num = len(bub_dict["bub_num"])

    bub_num = start_num
    bubs_found = False
    first_found = False
    stop_frame = 100

    #get number of frames in event
    keys = ev['cam'][f'c{cam}'].keys()
    frames = 0
    for key in keys:
        if 'frame' in key:
            frames+=1

    for i in range(1,frames-1):

        if i>stop_frame:
            break
    
        im_num = i
    
        thisIm = np.float32(np.average(ev['cam'][f'c{cam}'][f'frame{im_num-1}'],axis=2))
        nextIm = np.float32(np.average(ev['cam'][f'c{cam}'][f'frame{im_num}'],axis=2)) #may have to be changed later to account for missing frames
    
        #get diff and mask out noise
        preMask_diff = abs(thisIm-nextIm)
        preMask_diff[preMask_diff<noise_thresh] = 0
        
        #apply mask to isolate bubble region depending on which camera we're viewing from
        diff = np.zeros((imShape[0],imShape[1]))
        diff[mask_circy,mask_circx] = preMask_diff[mask_circy,mask_circx]
        diff-=dip.GetSinglePixels(diff > 0)
    
        #find and label connected regions of nonzero pixels 
        #connectivity = 2 allows pixels to count as connected if they are diagonal from each other
        labeled = label(diff>0, connectivity = 2)
    
        #get properties of labeled regions
        props = regionprops(labeled, intensity_image = diff)
    
        #get 2.5 sigma threshold for region intensity
        intensities = [prop.intensity_mean for prop in props]
        intensity_thresh = np.average(intensities) + 2.5*np.std(intensities)
    
        #get areas so we can identify number of potential bubbles in image
        areas = np.array([prop.area for prop in props])
        area_mask = areas>=areas.max()/4
        largest_regions = np.array(props)[area_mask]
    
        #get region with most connected pixels 
        largest_region = max(largest_regions, key=lambda prop: prop.area)
    
        #see if largest region meets area and intensity thresholds
        #if yes, or if we saw a bubble in a previous frame, we have a bubble
        if (largest_region.area>=10 and largest_region.intensity_mean>=intensity_thresh) or bubs_found==True:
    
            if bubs_found==False:
                stop_frame = im_num+5
                bubs_found = True
                first_found = True
                
            #estimate radii
            min_est_rad = np.round((largest_region.axis_minor_length+largest_region.axis_major_length)/4)
            if min_est_rad - 3 <=3:
                min_rad = 3
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
                    this_layer += np.roll(diff, offset,(0,1))
    
                if rad==rad_cands[0]:
                    accum = this_layer
                else:
                    accum = np.dstack((accum, this_layer))
    
            accum_shape = accum.shape
    
            #get vote threshold -- 80% of highest peak number of votes
            largest_cy, largest_cx = largest_region.centroid
            largest_y, largest_x = disk((largest_cy, largest_cx), 20, shape=imShape)
            max_votes = np.max(accum[largest_y, largest_x])
            vote_thresh = max_votes*.8
            
            #now we enter the candidate loop and decide based on votes whether to keep or discard each candidate
            for region in largest_regions:
           
                #helpful check for frames with no bubbles
                if region.axis_minor_length== 0 or region.axis_major_length==0 or region.area<10:
                    continue
                
                #if this is the first frame in which we see a bubble, scan back some frames to get true t0
                if first_found==True:
                    
                    #look for bubble in relevant region in past images
                    t0_found = False
                    j = 1
                    while t0_found==False:
                        
                        #get diff
                        pastIm = np.float32(np.average(ev['cam'][f'c{cam}'][f'frame{im_num-j}'], axis=2))
                        pastDiff = abs(pastIm - refIm)
                        pastDiff[pastDiff<noise_thresh] = 0
                        pastDiff-=dip.GetSinglePixels(pastDiff > 0)
    
                        #perform CHT
                        rad_cands = [2,3,4]
                        for rad in rad_cands:
                            circx, circy = circle_perimeter(600, 400, int(rad), shape=Tshape)
                            dx = circx-600
                            dy = circy-400
                            offsets = [(dx[i],dy[i]) for i in range(len(dx))]
    
                            this_layer = np.zeros((imShape[0], imShape[1]))
                            for offset in offsets:
                                this_layer += np.roll(pastDiff, offset,(0,1))
    
                            if rad==rad_cands[0]:
                                past_accum = this_layer
                            else:
                                past_accum = np.dstack((accum, this_layer))
    
                        past_accum_shape = past_accum.shape
    
                        #paste relevant region of accumulator onto blank image
                        pastregIm = np.zeros(past_accum_shape)
                        pastregIm[largest_y, largest_x] = past_accum[largest_y, largest_x]
    
                        #find peak candidate 
                        pcy, pcx, rad_ind = np.unravel_index(np.argmax(pastregIm), past_accum_shape)
                        prad = rad_ind + rad_cands[0]
    
                        #check if bubble candidate meets intensity thresh for the noise in the image
                        past_intensity_thresh = np.average(pastDiff) + 2.5*np.std(pastDiff)
                        bub_coords_y, bub_coords_x = disk((pcy, pcx), prad, shape=imShape)
                        if np.average(pastDiff[bub_coords_y, bub_coords_x]) >= past_intensity_thresh:
                            #add bubble to dictionary
                            bub_dict["bub_num"].append([bub_num])
                            bub_dict["cam"].append([cam])
                            bub_dict["pos"].append([pcx, pcy])
                            bub_dict["radius"].append([prad])
                            bub_dict["significance"].append([1.0])
                            bub_dict["frame"].append([im_num-j])

                            bub_num+=1
                            j+=1
                    
                        else:
                            #if no bubble, don't bother scanning through remaining frames
                            t0_found = True
                    
                    first_found = False
                
                #now continue on finding bubbles in successive frames as usual
                #get number of votes in this region of the accumulator
                rcy, rcx = region.centroid
                ry, rx = disk((rcy,rcx), 20, shape=imShape)
                regIm = np.zeros(accum_shape)
                regIm[ry,rx] = accum[ry, rx]
                votes = np.max(regIm)
                
                if votes>=vote_thresh:
                    
                    #get peak candidate in this constrained region of connected pixels
                    pcy, pcx, rad_ind = np.unravel_index(np.argmax(regIm), accum_shape)
                    prad = rad_ind + min_rad
    
                    #add bubble to dictionary
                    bub_dict["bub_num"].append([bub_num])
                    bub_dict["cam"].append([cam])
                    bub_dict["pos"].append([pcx, pcy])
                    bub_dict["radius"].append([prad])
                    bub_dict["significance"].append([votes/max_votes])
                    bub_dict["frame"].append([im_num-j])

                    bub_num+=1
                
                    #zero out accumulator array for all radii around this point
                    circy, circx = disk((pcy, pcx), prad+num_pix_in_neighborhood, shape=imShape)
                    accum[circy,circx] = 0

    return bub_dict

def BubbleFinder(ev, num_pix_in_neighborhood = 20, noise_thresh = 10):
    
    out = _new_bub_dict()
    out = FindBubbles(ev, 1, num_pix_in_neighborhood, noise_thresh, bub_dict=out)
    out = FindBubbles(ev, 2, num_pix_in_neighborhood, noise_thresh, bub_dict=out)
    out = FindBubbles(ev, 3, num_pix_in_neighborhood, noise_thresh, bub_dict=out)

    if len(out["bub_num"]) == 0:
        raise ValueError("No bubbles found in event")

    return out
