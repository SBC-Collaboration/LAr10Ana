import numpy as np
from scipy import ndimage
from skimage.draw import circle_perimeter
from skimage.draw import disk
import cv2 as cv

#output dictionary structure:
#camera --> frame --> bubble --> individual bubble RQs (position, radius, goodness of fit)
#for example, dict['cam3']['frame7']['bubble2']['radius'] gives radius of 3rd bubble found in 7th frame taken by camera 3 in the input event

#inputs for GetBubbles:
#ev: event loaded from GetEvent
#cam_list: list of cameras active in this event
#noise_thresh = 15: values below this threshold in the diff will be masked out
#num_pix_in_neighborhood = 10: radius in pixels around a prominent peak in the accumulator array set to zero before searching for next most prominent peak
#gof_thresh = 20: threshold below which goodness of fit metric must fall to declare a bubble
#all defaults very rough at this point, will tune further in the future

def GetGoF(white_x, white_y, center_y, center_x, radius, shape):
    
    #get pixels in the circles found by the CHT
    circy, circx = circle_perimeter(center_y, center_x, radius, shape=shape)
    #gof = goodness of fit
    gof_sum = 0
    num_pix = len(circy)
    for i in range(num_pix):
        #find distances between each white pixel and this pixel
        distances = np.sqrt((white_y - circy[i])**2 + (white_x - circx[i])**2)
        gof_sum+=distances.min()**2 #keep square of the distance between this pixel and the nearest white pixel
        
    return gof_sum/num_pix/radius #normalize these scores by radius so large asymmetric 
                                  #bubbles don't have dramatically larger scores than tiny bubbles

def FindBubbles(ev,cam,noise_thresh = 15,num_pix_in_neighborhood = 10,gof_thresh = 20):

    keys = ev['cam'][f'c{cam}'].keys()
    frames = 0
    for key in keys:
        if 'frame' in key:
            frames+=1
    
    cam_dict = {}
    
    bubs_found = False

    for i in range(frames-1):

        im_num = i+1
    
        #get image
        thisIm = np.float32(np.average(ev['cam'][f'c{cam}'][f'frame{i}'],axis=2))
        nextIm = np.float32(np.average(ev['cam'][f'c{cam}'][f'frame{im_num}'],axis=2))
        imShape = nextIm.shape

        #replace this part of the code when I can more robustly differentiate noise/reflections from small bubbles
        '''
        if bubs_found==False:
            #get diff, mask out noise
            diff = abs(nextIm - thisIm)
            diff[(diff<=noise_thresh) | (diff>=255-noise_thresh)] = 0
        else:
            blankNextIm = np.zeros(imShape)
            blankThisIm = np.zeros(imShape)
            #only worry about the neighborhood around previously found bubbles
            for center_y, center_x, radius in zip(prev_cy, prev_cx, prev_rad):
                circy, circx = disk((np.int32(center_y), np.int32(center_x)), np.int32(radius)+10, shape=imShape)
                blankNextIm[circy,circx] = nextIm[circy,circx]
                blankThisIm[circy,circx] = thisIm[circy,circx]
            #now get diff, mask out noise
            diff = abs(blankNextIm - blankThisIm)
            diff[(diff<=noise_thresh) | (diff>=255-noise_thresh)] = 0
        '''
        diff = abs(nextIm - thisIm)
        diff[(diff<=noise_thresh) | (diff>=255-noise_thresh)] = 0
        
        #get rid of pixels with fewer than 3 nonzero neighbors
        nonzero = np.array(np.where(diff>0))
        diff_copy = np.copy(diff)

        for i in range(len(nonzero[0])):
            pixx, pixy = nonzero[0][i], nonzero[1][i]

            try:
                surrounding_pixels = 0
    
                if diff_copy[pixx+1,pixy]>0:
                    surrounding_pixels+=1
                if diff_copy[pixx,pixy+1]>0:
                    surrounding_pixels+=1
                if diff_copy[pixx,pixy-1]>0:
                    surrounding_pixels+=1
                if diff_copy[pixx-1,pixy]>0:
                    surrounding_pixels+=1
                if diff_copy[pixx+1,pixy+1]>0:
                    surrounding_pixels+=1
                if diff_copy[pixx+1,pixy-1]>0:
                    surrounding_pixels+=1
                if diff_copy[pixx-1,pixy-1]>0:
                    surrounding_pixels+=1
                if diff_copy[pixx-1,pixy+1]>0:
                    surrounding_pixels+=1
        
                if surrounding_pixels<3:
                    diff[pixx,pixy] = 0
                    
            except:
                diff[pixx,pixy] = 0 #accounts for when we roll off the image; we won't see bubbles at these edges anyhow
        
        #get radii candidates and prepare accumulator array
        
        if bubs_found==False:
            min_rad_cand = 3
            max_rad_cand = 10
        else:
            min_rad_cand = np.min(prev_rad)+1
            max_rad_cand = np.max(prev_rad)+10
        radii_cands = np.int32(np.arange(min_rad_cand,max_rad_cand,1))
        accum_shape = (imShape[0],imShape[1],len(radii_cands))
        

        #perform the CHT
        accum = np.zeros(accum_shape)
        edge_pix = np.where(diff!=0)
        x = edge_pix[1]
        y = edge_pix[0]
        Tshape = thisIm.T.shape
        for i in range(x.shape[0]):
            for k in range(len(radii_cands)):
                xi = x[i]
                yi = y[i]
                circx, circy = circle_perimeter(xi, yi, radii_cands[k], shape=Tshape)
                accum[circy,circx,k] += diff[yi][xi]
                
        #identify peak candidates
        peaks = np.where(accum>0) 
        cy = peaks[0]
        cx = peaks[1]
        radii = radii_cands[peaks[2]]

        if bubs_found==False:
            max_bub = 5
        else:
            max_bub = len(prev_rad)

        #get peaks with highest accumulator score
        circle_thresh = np.max(accum)/2 #not accepting low quality circle cands

        bub_cx = []
        bub_cy = []
        bub_rad = []

        for i in range(max_bub):
            if np.max(accum)>circle_thresh:
                this_peak = np.where(accum==accum.max())
                this_cy = this_peak[0][0]
                this_cx = this_peak[1][0]
                this_rad = radii_cands[this_peak[2][0]]

                bub_cx.append(this_cx)
                bub_cy.append(this_cy)
                bub_rad.append(this_rad)

                circy, circx = disk((this_cy, this_cx), this_rad+num_pix_in_neighborhood, shape=imShape)
                accum[circy,circx] = 0

        #get goodness of fit for each circle
        diff[diff>0] = 255
        white = np.where(diff==255)
        white_y = white[0]
        white_x = white[1]
        gofs = []
        for i in range(len(bub_rad)):
            gof = GetGoF(white_x,white_y,bub_cy[i],bub_cx[i],bub_rad[i],imShape)
            gofs.append(gof)
        gofs = np.array(gofs)

        #mask out circles with bad gof
        gof_mask = gofs<gof_thresh
        gofs = gofs[gof_mask]
        bub_rad = np.array(bub_rad)
        bub_cx = np.array(bub_cx)
        bub_cy = np.array(bub_cy)
        
        bub_rad = bub_rad[gof_mask]
        bub_cx = bub_cx[gof_mask]
        bub_cy = bub_cy[gof_mask]

        num_bubs_found = len(bub_rad)

        #save bubble info if there's info to save
        if num_bubs_found>0:
            bubs_found = True              
            frame_dict = {}
            for i in range(num_bubs_found):
                bub_dict = {}
                bub_dict['Pos'] = (np.uint32(bub_cx[i]),np.uint32(bub_cy[i]))
                bub_dict['Radius'] = bub_rad[i]
                bub_dict['GOF'] = gofs[i]

                frame_dict[f'Bubble{i}'] = bub_dict
            cam_dict[f'Frame{im_num}'] = frame_dict
        
            #update circle params for next iteration and return dict
            prev_rad, prev_cx, prev_cy = bub_rad, bub_cx, bub_cy
        else:
            if im_num<15:
                #im_num>15 is a bandaid I'll have to revisit later after I can better distinguish
                #between noise/reflections and small bubbles
                bubs_found = False

    return cam_dict

def GetBubbles(ev,cam_list,noise_thresh = 15,num_pix_in_neighborhood = 10,gof_thresh = 20):

    ev_dict = {}
    if 1 in cam_list:
        cam1_dict = FindBubbles(ev,1,noise_thresh,num_pix_in_neighborhood,gof_thresh)
        ev_dict['Cam1'] = cam1_dict
    if 2 in cam_list:
        cam2_dict = FindBubbles(ev,2,noise_thresh,num_pix_in_neighborhood,gof_thresh)
        ev_dict['Cam2'] = cam2_dict
    if 3 in cam_list:
        cam3_dict = FindBubbles(ev,3,noise_thresh,num_pix_in_neighborhood,gof_thresh)
        ev_dict['Cam3'] = cam3_dict

    return ev_dict
