from sbcbinaryformat import Streamer, Writer
import numpy as np
from sklearn.cluster import DBSCAN

def CleanBubDict(bub_dict):

    #keep dictionary entries starting at 0 and going in order even if we skip a bubble
    nGoodBubs = 0

    #set up dictionary entry for each bubble being tracked
    out = {}
    
    #for testing purposes
    #run_nbub = []

    for ev in range(0,bub_dict['ev'].max()+1):
        
        #get event info
        evt = bub_dict['ev'] == ev
        cams = bub_dict['cam'][evt]
        frames = bub_dict['frame'][evt]
        pos = bub_dict['pos'][evt]
        rads = bub_dict['radius'][evt]
        sigs = bub_dict['significance'][evt]
        run = bub_dict['runid'][evt]

        #for testing purposes
        #ev_nbub = []
        
        for cam in range(1,4):

            #get info for this cam
            this_cam = cams == cam
            cam_n = cams[this_cam]
            frames_n = frames[this_cam]
            pos_n = pos[this_cam]
            rads_n = rads[this_cam]
            sigs_n = sigs[this_cam]
            run_n = run[this_cam]

            if len(frames_n)==0:
                continue

            #use dbscan to estimate number of clusters/ bubbles
            db = DBSCAN(eps=10, min_samples=2).fit(pos_n)
            db_labels = db.labels_
            n_db_clusters = len(set(db_labels)) - (1 if -1 in db_labels else 0)
            best_clusters = db_labels
            nbub = n_db_clusters
            
            #if a cluster has less spread in its radii than the radii for all detections in this cam,
            #it's more likely one bubble has been broken into multiple clusters; we expect each bubble to 
            #grow during an event and therefore have a range of radii; group this cluster with the nearest one,
            #reassign cluster assignments as necessary, and adjust nbub estimate
            if nbub>1:
                grouped = 0
                cents = np.array([(np.average(pos_n[best_clusters==i][:,0]), np.average(pos_n[best_clusters==i][:,1])) for i in range(nbub)])
                centsx, centsy = cents[:,0], cents[:,1]
                distsx = abs(np.subtract.outer(centsx, centsx))**2
                distsy = abs(np.subtract.outer(centsy, centsy))**2
                dists = np.sqrt(distsx+distsy)
                dists[dists==0] = 100
                for i in range(nbub):
                    if np.std(rads_n[best_clusters==i])<np.std(rads_n) and nbub>=1: 
                        best_clusters[best_clusters==i] = np.argmin(dists[i])
                        dists[i,:] = 100
                        dists[:,i] = 100
                        grouped+=1
                nbub-=grouped
    
                #reset best_clusters to start at 0, as rest of code requires this
                clean_nbub = 0
                for val in np.unique(best_clusters):
                    if val==-1:
                        continue
                    best_clusters[best_clusters==val] = clean_nbub
                    clean_nbub+=1
    
            #identify bubbles that don't occur in 2+ frames with others (or at all), AKA false noise tags
            all_frames = [frames_n[best_clusters==i] for i in range(nbub)]
            intersections = np.zeros((nbub,nbub))
            for i in range(nbub):
                for j in range(nbub):
                    if i==j:
                        intersections[i][j] = len(all_frames[i])
                    else:
                        intersections[i][j] = len(np.intersect1d(all_frames[i],all_frames[j]))
            vals, counts = np.unique(np.where(intersections<2)[0], return_counts=True)
            bad_bubs = vals[counts>1]

            #for testing only
            #ev_nbub.append(nbub - len(badbubs))

            #set up dictionary entry for each bubble being tracked
            cam_out = {f'Bub{i}': {} for i in range(nGoodBubs, nGoodBubs + nbub - len(bad_bubs))}
            #cam_out = {}
    
            #save grouped bubble info
            for i in range(nbub):
            
                #don't add if this is a talse tag we've identified
                if i in bad_bubs:
                    continue
            
                #get frames and positions associated with this bubble
                f = frames_n[best_clusters==i]
                p = pos_n[best_clusters==i]
                r = rads_n[best_clusters==i]
                s = sigs_n[best_clusters==i]
                ri = run_n[best_clusters==i]
            
                #remove double tags within the same frame; keep whichever tag is closest to cluster centroid
                centroid = (np.average(p[:,0]), np.average(p[:,1]))
                f_clean = np.unique(f)
                p_clean = []
                r_clean = []
                s_clean = []
                ri_clean = []
                for frame in f_clean:
                    frame_mask = f==frame
                    tags = p[frame_mask]
                    tag_rads = r[frame_mask]
                    tag_sigs = s[frame_mask]
                    tag_rid = ri[frame_mask]
                    if tags.size==1:
                        p_clean.append(tags)
                        r_clean.append(tag_rads)
                        s_clean.append(tag_sigs)
                        ri_clean.append(tag_rid)
                    else:
                        dists = [np.sqrt((tag[0]-centroid[0])**2 + (tag[1]-centroid[1])**2) for tag in tags]
                        indx = np.argmin(dists)
                        p_clean.append(tags[indx])
                        r_clean.append(tag_rads[indx])
                        s_clean.append(tag_sigs[indx])
                        ri_clean.append(tag_rid[indx])

                
                cam_out[f'Bub{nGoodBubs}']['frames'] = f_clean
                cam_out[f'Bub{nGoodBubs}']['pos'] = p_clean
                cam_out[f'Bub{nGoodBubs}']['rad'] = r_clean
                cam_out[f'Bub{nGoodBubs}']['cam'] = np.full(len(f_clean),cam)
                cam_out[f'Bub{nGoodBubs}']['ev'] = np.full(len(f_clean),ev)
                cam_out[f'Bub{nGoodBubs}']['significance'] = s_clean
                cam_out[f'Bub{nGoodBubs}']['runid'] = ri_clean
                            
                nGoodBubs+=1

            out.update(cam_out)
        #run_nbub.append([ev,max(ev_nbub)])

    return out #, run_nbub