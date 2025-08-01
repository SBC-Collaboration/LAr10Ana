import numpy as np

SAMPLE_FREQ = 62.5 # MHz

def nSampAvg(wvf,n):
    arr_list = []
    for i in range(n):
        if i==n-1:
            arr_list.append(wvf[i:][:][:])
        else:
            arr_list.append(wvf[i:-(n-1-i)][:][:])
    stack_arr = np.stack(arr_list,axis=-1)
    return np.average(stack_arr,axis=-1)

def FindEdges(traces,peak_inds,thresh,N):
    #from peak, work backwards in time until signal has remained under the threshold for N samples
    #simultaneously, work forwards in time until signal dips back below threshold

    t0 = np.zeros(traces.shape[1]) #pre-allocating memory instead of appending, is this the right idea?
    tf = np.zeros(traces.shape[1])
    for i in range(traces.shape[1]):
        peak_ind = peak_inds[i]
        #if no pulse at this index, leave t0 and tf as zero and skip to next index
        if peak_ind==0:
            continue
        #scanning backwards
        j = 1
        samps_under_thresh = 0
        while samps_under_thresh<N:
            if peak_ind-j<=2:
                j = peak_ind
                break 
            if traces[peak_ind-j,i]<=thresh[i]:
                samps_under_thresh+=1
            else:
                samps_under_thresh = 0
            j+=1
        t0[i] = peak_ind-j
        #scanning forwards
        k = 1
        samps_under_thresh = 0
        while samps_under_thresh<N: #can change this later with nicer data if the tails are long
            if peak_ind+k>traces.shape[0]-2:
                k = traces.shape[0]-peak_ind
                break
            if traces[peak_ind+k,i]<=thresh[i]:
                samps_under_thresh+=1
            else:
                samps_under_thresh = 0
            k+=1
        tf[i] = peak_ind+k
    return t0, tf

def aftXX(traces,t0,areas,XX): #XX given as a percentage
    aft = np.zeros(traces.shape[1])
    for i in range(traces.shape[1]):
        trace_cumsum = np.cumsum(traces[:,i])
        target_area = areas[i]*XX*0.01
        diffs = trace_cumsum - target_area
        aft[i] = np.argmin(abs(diffs))
    return aft #returns index at which hit area is XX% of total hit area

#arguments are: event, number of samples for which trace must remain under rms*5 to declare a pulse edge, number of samples used in n-sample averaging for #waveform smoothing, whether or not you want to use a threshold in order to declare that the highest peak in a trace is really a pulse (set to True when #making finger plots or otherwise looking at pulses on the order of a single photon)
def SiPMPulses(ev,N,samp,SiPM_cal=False):
    # Stuff to save, with defaults
    default_output = dict(
        baseline=np.array([]),
        rms=np.array([]),
        hit_t0=np.array([]),
        hit_t100=np.array([]),
        hit_t10=np.array([]),
        hit_t50=np.array([]),
        hit_t90=np.array([]),
        hit_area_ch=np.array([]),
        hit_area_tot=np.array([]),
        hit_amp_ch=np.array([]),
        hit_amp_tot=np.array([]),
        second_pulse_t0=np.array([]),
        second_pulse_amp=np.array([]), #can easily add more info about second pulse if desired
        num_pulses=np.array([]),
        summed_waveforms=np.array([])
    )
    out = default_output

    if ev is None:
        return out

    # One event has N SiPMs readout M times, each with T samples
    # Each readout is stored is ev["sipm_traces"], with a shape (M, N, T)

    # For each SiPM, for each readout, obtain the t0 and the voltage

    # Waveform in (ticks, ADC)
    traces = ev["scintillation"]["Waveforms"].T.astype(float)

    # Subtract offset and convert to mV
    for i_sipm in range(traces.shape[1]):
        group = i_sipm // 8
        chan  = i_sipm % 8
        group_ctrl = ev["run_control"]["caen"]["group%i" % group]
        offset = group_ctrl["offset"] + group_ctrl["ch_offset"][chan]
        range_mV = float(group_ctrl["range"][:-4])

        traces[:, i_sipm, :] -= offset
        traces[:, i_sipm, :] *= range_mV


    decimation = ev["run_control"]["caen"]["global"]["decimation"]
    sample_rate =  SAMPLE_FREQ/(2**decimation)

    # obtain the leading baseline and RMS
    N_SAMPLE_BASELINE = 40 

    baseline = traces[:N_SAMPLE_BASELINE].mean(axis=0)
    rms = traces[:N_SAMPLE_BASELINE].std(axis=0)

    # flip the trace and correct for baseline
    raw_traces = -(traces - baseline)

    #smooth waveforms
    smooth_traces = nSampAvg(raw_traces,samp)

    #find peaks for first pulse (if there is one)
    N_SIGMA_THRESHOLD = 5
    rms_thresh = rms*N_SIGMA_THRESHOLD
    
    peaks = np.max(smooth_traces,axis=0)
    peak_inds = np.argmax(smooth_traces,axis=0)

    threshMask = peaks<rms_thresh
    if SiPM_cal==False:
        peaks[threshMask] = 0
        peak_inds[threshMask] = 0

    #now find pulse edges
    t0 = np.zeros((smooth_traces.shape[1],smooth_traces.shape[2]))
    tf = np.zeros((smooth_traces.shape[1],smooth_traces.shape[2]))
    for ch in range(smooth_traces.shape[1]):
        ch_traces = smooth_traces[:,ch,:]
        ch_traces[:,threshMask[ch]] = 0
        ch_t0, ch_tf = FindEdges(ch_traces,peak_inds[ch,:],rms_thresh[ch,:],10)
        t0[ch] = ch_t0
        tf[ch] = ch_tf

    #now check whether there's a second pulse in the trace
    for ch in range(smooth_traces.shape[1]):
        for i in range(smooth_traces.shape[2]):
            smooth_traces[np.int32(t0[ch,i]):np.int32(tf[ch,i]),ch,i] = 0

    peaks2 = np.max(smooth_traces,axis=0)
    peak_inds2 = np.argmax(smooth_traces,axis=0)

    threshMask2 = peaks2<rms_thresh
    peaks2[threshMask2] = 0
    peak_inds2[threshMask2] = 0

    t02 = np.zeros((smooth_traces.shape[1],smooth_traces.shape[2]))
    tf2 = np.zeros((smooth_traces.shape[1],smooth_traces.shape[2]))
    for ch in range(smooth_traces.shape[1]):
        ch_traces = smooth_traces[:,ch,:]
        ch_traces[:,threshMask2[ch]] = 0
        ch_t0, ch_tf = FindEdges(ch_traces,peak_inds2[ch,:],rms_thresh[ch,:],10)
        t02[ch] = ch_t0
        tf2[ch] = ch_tf 
    
    #record how many pulses found in each scintillation trigger in each channel
    peaklets = np.zeros((smooth_traces.shape[1],smooth_traces.shape[2]))
    for ch in range(smooth_traces.shape[1]):
        peaklets[ch][np.where(tf[ch]!=0)[0]] = 1
        peaklets[ch][np.where(tf2[ch]!=0)[0]] = 2 #if I were spending more time on this, I'd be cleverer about not updating some of these entries twice

    #delete traces outside of found pulses and sum 
    #also switching to raw traces now since these are the traces we'll use to calculate the rest of our RQs
    #also also getting areas by channel
    areas_ch = np.zeros((raw_traces.shape[1],raw_traces.shape[2]))
        #should be noted that I guess all areas are in mV*samples right now, can change x dimension to ns later if desired
        #and also of course we can change mV to phd once calibrations are done

    for ch in range(raw_traces.shape[1]):
        for i in range(raw_traces.shape[2]):
            if peaklets[ch,i]==0:
                raw_traces[:,ch,i] = 0
            elif peaklets[ch,i]==1:
                raw_traces[:np.int32(t0[ch,i]),ch,i] = 0
                raw_traces[np.int32(tf[ch,i]):,ch,i] = 0
            elif peaklets[ch,i]==2:
                if np.int32(t0[ch,i])<np.int32(t02[ch,i]):
                    raw_traces[:np.int32(t0[ch,i]),ch,i] = 0
                    raw_traces[np.int32(tf[ch,i]):np.int32(t02[ch,i]),ch,i] = 0
                    raw_traces[np.int32(tf2[ch,i]):,ch,i] = 0
                else:
                    raw_traces[:np.int32(t02[ch,i]),ch,i] = 0
                    raw_traces[np.int32(tf2[ch,i]):np.int32(t0[ch,i]),ch,i] = 0
                    raw_traces[np.int32(tf[ch,i]):,ch,i] = 0
                
            areas_ch[ch,i] = np.sum(raw_traces[:,ch,i])
        
    summed_traces = np.sum(raw_traces,axis=1)

    #Awesome!  Now let's get the rest of our RQs :]
    areas_tot = np.sum(summed_traces,axis=0)

    t0_tot = np.zeros(summed_traces.shape[1])
    for i in range(summed_traces.shape[1]):
        try:
            t0_tot[i] = np.where(summed_traces[:,i]!=0)[0][0]
        except:
            continue

    t100 = aftXX(summed_traces,t0_tot,areas_tot,100)
    t90 = aftXX(summed_traces,t0_tot,areas_tot,90)
    t10 = aftXX(summed_traces,t0_tot,areas_tot,10)
    t50 = aftXX(summed_traces,t0_tot,areas_tot,50)

    amps = np.max(summed_traces,axis=0)

    #finally, save our RQs
    out["baseline"] = baseline
    out["rms"] = rms
    out["hit_t0"] = t0_tot #/ sample_rate
    out["hit_t100"] = t100 #/ sample_rate
    out["hit_t10"] = t10 #/ sample_rate
    out["hit_t50"] = t50 #/ sample_rate
    out["hit_t90"] = t90 #/ sample_rate
    out["hit_area_ch"] = areas_ch
    out["hit_area_tot"] = areas_tot
    out["hit_amp_ch"] = peaks
    out["hit_amp_tot"] = amps
    out["second_pulse_t0"] = t02 #/ sample_rate
    out["second_pulse_amp"] = peaks2
    out["num_pulses"] = peaklets
    out["summed_waveforms"] = summed_traces

    #got rid of "[t0==0] = nan" bit because we can use num_pulses to mask if desired
    
    return out

if __name__ == "__main__":
    from GetEvent import GetEvent

    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250611_1/"
    TEST_EVENT = 0
    out = SiPMPulses(GetEvent(TEST_RUN, TEST_EVENT))
    print("OUTPUTS:", *out.keys())
    print("BASELINES:", ", ".join(map(str, list(out["baseline"].mean(axis=-1)))))
    print("RMS:", ", ".join(map(str, list(out["rms"].mean(axis=-1)))))
    print("NHIT:", ", ".join(map(str, list((~np.isnan(out["hit_t0"])).astype(int).sum(axis=-1)))))
