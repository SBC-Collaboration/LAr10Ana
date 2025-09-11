import numpy as np
from scipy import ndimage
from skimage.feature import canny
from skimage.draw import circle_perimeter
from skimage.draw import disk
import cv2 as cv

#output dictionary structure:
#camera --> frame --> bubble --> individual bubble RQs (position, radius, goodness of fit)
#for example, dict['cam3']['frame7']['bubble2']['radius'] gives radius of 3rd bubble found in 7th frame taken by camera 3 in the input event

def GetGoF(white_x, white_y, center_y, center_x, radius, shape):
    
    #get pixels in the circles found by the CHT
    circy, circx = circle_perimeter(center_y, center_x, radius, shape=shape)
    #gof = goodness of fit
    gof_sum = 0
    num_pix = len(circy)
    for i in range(num_pix):
        #find distances between each white pixel and this pixel
        distances = np.sqrt((white_y - circy[i])**2 + (white_x - circx[i])**2)
        gof_sum+=distances.min()**2
        
    return gof_sum/num_pix

def FindBubbles(ev,cam,noise_thresh = 30,rad_thresh = 15,sig_high = 7,sig_low = 3,num_pix_in_neighborhood = 5,gof_thresh = 110):

    keys = ev['cam'][f'c{cam}'].keys()
    frames = 0
    for key in keys:
        if 'frame' in key:
            frames+=1

    cam_dict = {}
    
    bubs_found = False

    #get reference image
    thisIm = np.float32(np.average(ev['cam'][f'c{cam}']['frame0'],axis=2))

    for i in range(frames-1):

        im_num = i+1
    
        #get image
        nextIm = np.float32(np.average(ev['cam'][f'c{cam}'][f'frame{im_num}'],axis=2))
        imShape = nextIm.shape
        
        if bubs_found==False:
            #get diff, mask out noise, turn into binary image, apply Canny edge detection
            diff = abs(nextIm - thisIm)
            diff[(diff<=noise_thresh) | (diff>=255-noise_thresh)] = 0
            diff[diff>0] = 255
            edges = canny(diff, sigma=sig_low, low_threshold=10, high_threshold=50)
        else:
            blankNextIm = np.zeros(imShape)
            #only worry about the neighborhood around previously found bubbles
            for center_y, center_x, radius in zip(prev_cy, prev_cx, prev_rad):
                circy, circx = disk((np.int32(center_y), np.int32(center_x)), np.int32(radius)+10, shape=imShape)
                blankNextIm[circy,circx] = nextIm[circy,circx]
            #now get diff, mask out noise, turn into binary image, apply Canny edge detection
            diff = abs(nextIm - thisIm)
            diff[(diff<=noise_thresh) | (diff>=255-noise_thresh)] = 0
            diff[diff>0] = 255
            if np.max(prev_rad)>rad_thresh: 
                edges = canny(diff, sigma=sig_high, low_threshold=10, high_threshold=50) #higher sigma better for larger, less round bubbles
            else:
                edges = canny(diff, sigma=sig_low, low_threshold=10, high_threshold=50)
                
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
        edge_pix = np.where(edges==True)
        x = edge_pix[1]
        y = edge_pix[0]
        Tshape = edges.T.shape
        for i in range(x.shape[0]):
            for k in range(len(radii_cands)):
                circx, circy = circle_perimeter(x[i], y[i], radii_cands[k], shape=Tshape)
                accum[circy,circx,k] += 1
                
        #identify peaks
        peaks = np.where(accum>0)
        cy = peaks[0]
        cx = peaks[1]
        radii = radii_cands[peaks[2]]

        if bubs_found==False:
            max_bub = 5
        else:
            max_bub = len(prev_rad)
            
        circle_thresh = np.max(accum)/3 #not accepting low quality circle cands

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
        white = np.where(np.uint8(edges)==1)
        white_y = white[0]
        white_x = white[1]
        gofs = []
        for i in range(len(bub_rad)):
            gof = GetGoF(white_x,white_y,bub_cy[i],bub_cx[i],bub_rad[i],imShape)
            gofs.append(gof)
        gofs = np.array(gofs)

        #figure out whether we've actually got bubbles
        if bubs_found==False: #can't use on later images/ larger bubbles because they're less round and thus have worse fits
            gof_mask = gofs<gof_thresh
            gofs = gofs[gof_mask]

            bub_rad = np.array(bub_rad)
            bub_cx = np.array(bub_cx)
            bub_cy = np.array(bub_cy)
        
            bub_rad = bub_rad[gof_mask]
            bub_cx = bub_cx[gof_mask]
            bub_cy = bub_cy[gof_mask]

        num_bubs_found = len(bub_rad)
        if num_bubs_found>0:
            bubs_found = True
            #save bubble info
            frame_dict = {}
            for i in range(num_bubs_found):
                bub_dict = {}
                bub_dict['Pos'] = (np.uint32(bub_cx[i]),np.uint32(bub_cy[i]))
                bub_dict['Radius'] = bub_rad[i]
                bub_dict['GOF'] = gofs[i]

                frame_dict[f'Bubble{i}'] = bub_dict
            cam_dict[f'Frame{im_num}'] = frame_dict
        
        #update circle params for next iteration and return dict
        if bubs_found==True:
            prev_rad, prev_cx, prev_cy = bub_rad, bub_cx, bub_cy

    return cam_dict

def GetBubbles(ev,noise_thresh = 30,rad_thresh = 15,sig_high = 7,sig_low = 3,num_pix_in_neighborhood = 5,gof_thresh = 110):

    ev_dict = {}

    cam1_dict = FindBubbles(ev,1,noise_thresh,rad_thresh,sig_high,sig_low,num_pix_in_neighborhood,gof_thresh)
    cam2_dict = FindBubbles(ev,2,noise_thresh,rad_thresh,sig_high,sig_low,num_pix_in_neighborhood,gof_thresh)

    ev_dict['Cam1'] = cam1_dict
    ev_dict['Cam2'] = cam2_dict

    return ev_dict