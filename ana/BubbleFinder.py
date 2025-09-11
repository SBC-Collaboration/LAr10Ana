import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from skimage.feature import canny
from skimage.draw import circle_perimeter
from skimage.draw import disk
from PIL import Image
import cv2 as cv
import os as os

def GetGoF(white_x, white_y, center_y, center_x, radius, shape):
    
    #get pixels in the circles found by the CHT
    circy, circx = circle_perimeter(center_y, center_x, radius, shape=shape)
    #gof = goodness of fit
    gof_sum = 0
    for i in range(len(circy)):
        #find distances between each white pixel and this pixel
        distances = np.sqrt((white_y - circy[i])**2 + (white_x - circx[i])**2)
        gof_sum+=distances.min()**2
        
    return gof_sum/len(circy)

def FindBubbles(camFileList, path, plotMode = False):

    cam_dict = {}
    
    bubs_found = False

    #get reference image
    thisImPath = path + camFileList[0]
    thisIm = np.array(Image.open(thisImPath))

    for i in range(len(camFileList)-1):

        im_num = i+1
    
        #get image
        nextImPath = path + camFileList[im_num]
        nextIm = np.array(Image.open(nextImPath))
        
        if bubs_found==False:
            #get diff, mask out noise, turn into binary image, apply Canny edge detection
            diff = nextIm - thisIm
            noise_thresh = 30
            diff[diff<=noise_thresh] = 0
            diff[diff>=255-noise_thresh] = 0
            diff[diff>0] = 255
            edges = canny(diff, sigma=3, low_threshold=10, high_threshold=50)
        else:
            blankNextIm = np.zeros(nextIm.shape)
            #only worry about the neighborhood around previously found bubbles
            for center_y, center_x, radius in zip(prev_cy, prev_cx, prev_rad):
                circy, circx = disk((np.int32(center_y), np.int32(center_x)), np.int32(radius)+10, shape=nextIm.shape)
                blankNextIm[circy,circx] = nextIm[circy,circx]
            #now get diff, mask out noise, turn into binary image, apply Canny edge detection
            diff = abs(nextIm - thisIm)
            diff = nextIm - thisIm
            noise_thresh = 30
            diff[diff<=noise_thresh] = 0
            diff[diff>=255-noise_thresh] = 0
            diff[diff>0] = 255
            if np.max(prev_rad)>15: #thresh to be fine tuned in future
                edges = canny(diff, sigma=7, low_threshold=10, high_threshold=50) #higher sigma better for larger, less round bubbles
            else:
                edges = canny(diff, sigma=3, low_threshold=10, high_threshold=50)
                
        #get radii candidates and prepare accumulator array
        if bubs_found==False:
            min_rad_cand = 3
            max_rad_cand = 10
        else:
            min_rad_cand = np.min(prev_rad)+1
            max_rad_cand = np.max(prev_rad)+10
        radii_cands = np.int32(np.arange(min_rad_cand,max_rad_cand,1))
        accum_shape = (edges.shape[0],edges.shape[1],len(radii_cands))

        #perform the CHT
        accum = np.zeros(accum_shape)
        edge_pix = np.where(edges==True)
        x = edge_pix[1]
        y = edge_pix[0]
        for i in range(x.shape[0]):
            for k in range(len(radii_cands)):
                circx, circy = circle_perimeter(x[i], y[i], radii_cands[k], shape=edges.T.shape)
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
        num_pix_in_neighborhood = 5

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

                circy, circx = disk((this_cy, this_cx), this_rad+num_pix_in_neighborhood, shape=edges.shape)
                accum[circy,circx] = 0

        #get goodness of fit for each circle
        gof_thresh = 110 #somewhat arbitrary for now
        white = np.where(edges.astype(np.uint8)==1)
        white_y = white[0]
        white_x = white[1]
        gofs = []
        for i in range(len(bub_rad)):
            gof = GetGoF(white_x,white_y,bub_cy[i],bub_cx[i],bub_rad[i],nextIm.shape)
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

        if len(bub_rad)>0:
            bubs_found = True
            #save bubble info
            frame_dict = {}
            for i in range(len(bub_rad)):
                bub_dict = {}
                bub_dict['Pos'] = (np.uint32(bub_cx[i]),np.uint32(bub_cy[i]))
                bub_dict['Radius'] = bub_rad[i]
                bub_dict['GOF'] = gofs[i]
                
                frame_dict[f'Bubble{i}'] = bub_dict
            cam_dict[f'Frame{im_num}'] = frame_dict
        
        #plot, if desired
        if plotMode==True:
            
            fig, ax = plt.subplots(1,2,figsize=(8,6))

            ax[1].imshow(nextIm,cmap='gray')
            ax[1].axis('off')

            if bubs_found==True:
                for center_y, center_x, radius in zip(bub_cy, bub_cx, bub_rad):
                    for i in range(5): #drawing circles 5 pixels wide for visibility
                        circy, circx = circle_perimeter(np.int32(center_y), np.int32(center_x), np.int32(radius+i), shape=nextIm.shape)
                        nextIm[circy,circx] = 255 
                ax[0].scatter(bub_cx,bub_cy,c='r',s=1)
 
            ax[0].imshow(nextIm,cmap='gray')
            ax[0].axis('off')
            plt.show()
            
        #update circle params for next iteration and return dict
        if bubs_found==True:
            prev_rad, prev_cx, prev_cy = bub_rad, bub_cx, bub_cy

    return cam_dict

def GetImages(path):
    
    filesInPath = list(os.walk(path))[0][2]
    cam0FileList = []
    cam1FileList = []
    for file in filesInPath:
        if 'cam0image' in file:
            cam0FileList.append(file)
        elif 'cam1image' in file:
            cam1FileList.append(file)
    cam0FileList = np.sort(cam0FileList)
    cam1FileList = np.sort(cam1FileList)

    return cam0FileList, cam1FileList

def GetBubbles(path, plotMode = False):

    ev_dict = {}
    
    cam0, cam1 = GetImages(path)

    cam0_dict = FindBubbles(cam0, path, plotMode)
    cam1_dict = FindBubbles(cam1, path, plotMode)

    ev_dict['Cam0'] = cam0_dict
    ev_dict['Cam1'] = cam1_dict

    return ev_dict
