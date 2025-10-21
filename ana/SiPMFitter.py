# Queen's SiPM pulse fitting 

import sys

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 300
import scipy.optimize # For curve fitting

from scipy.optimize import curve_fit
from scipy import stats
import json

from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
from scipy.stats import norm

from ana  import BatchSiPMs

SAMPLE_FREQ = 62.5 # MHz

# Base functions

#This function just removes the droop, as written by Ben. Runs on a single trace
def droopFix(data, droop_tau=150,  # still just a guess
      t0=80):
    droop_A = droop_tau * (1-np.exp(-1./droop_tau))
    wf_corr = np.zeros(np.size(data))
    S = np.zeros(np.size(data)) # recursive correction factor
    for i in range(np.size(data)):
        if i < t0: # handles pre-trigger
            S[i] = 0
            wf_corr[i] = data[i]
        elif i == t0: # acquisition t0
            S[i] = 0
            wf_corr[i] = data[i]/droop_A
        else:
            S[i] = wf_corr[i-1] + np.exp(-1/droop_tau)*S[i-1]
            wf_corr[i] = data[i]/droop_A + (droop_A/droop_tau)*S[i]
    return(wf_corr)

#This is the denoising written by Gary Sweeney, runs on a single trace
def fft_denoise(signal, dt, cutoff_freq):
    # Get the fft and frequency of the waveform
    N = len(signal)
    fft_vals = np.fft.fft(signal)
    freqs = np.fft.fftfreq(N, dt)
    # Create a low-pass filter mask
    mask = np.abs(freqs) <= cutoff_freq
    filtered_fft = fft_vals * mask
    # Inverse FFT to get back to time domain
    denoised_signal = np.fft.ifft(filtered_fft).real
    return denoised_signal, freqs, np.abs(fft_vals), np.abs(filtered_fft)

#This is the function being used now, which seems to work well
def newFitFunc(t,*p): #time, and fit parameters (length that is dependent on number of pulses)
    # length of p (e.g. initial_guess) will tell us how many pulses to fit
    n = int(np.size(p)/5) #number of pulses to fit
    splitInput = np.array_split(np.array(p),n)  #splits the input p into n arrays of the components for each fit
    # nEXO SiPM pulse template, but with k = 0. arXiv:1903.03663
    # p[0] -> t0, p[1] -> pulse area (eventually proportional to nPE), p[2] -> baseline
    # p[3] -> tau_s (fall time), p[4] -> tau_r (rise time)
    f = np.heaviside((t-splitInput[0][0]),0.5)*splitInput[0][1]*(1/splitInput[0][3]*(np.exp(-(t-splitInput[0][0])/(splitInput[0][3]+splitInput[0][4])) - np.exp(-(t-splitInput[0][0])/(splitInput[0][4]))))
    if n>0:
        for i in range(1,n):
            f += np.heaviside((t-splitInput[i][0]),0.5)*splitInput[i][1]*(1/splitInput[i][3]*(np.exp(-(t-splitInput[i][0])/(splitInput[i][3]+splitInput[i][4])) - np.exp(-(t-splitInput[i][0])/(splitInput[i][4]))))
    return f  # f is now a sum of n calls of the single-pulse function

# L2 functions

##########
################
#This is the collection of functions that takes in the waveforms
#and outputs the transfomred waveforms and error
def prepareWaveforms(wfs, t0):
    range_V = 2 # from json
    dt = 16 # 1 sample = 16 ns (CAEN samples at 62.5 MHz), from json?
    cutoff = 0.005 # Cutoff at 200 MHz
    numEvents = np.shape(wfs)[0]
    numSipms = np.shape(wfs)[1]
    numSamples = np.shape(wfs)[2]
    #print('prepareWaveforms numEvents:',numEvents)
    #print('prepareWaveforms numSipms:',numSipms)
    wfFinal = []
    wfErrFinal = []
    #Start something to give uncertainty on these values
    for ev in range(numEvents):
        wfEvent = []
        wfErrEvent = []
        for sipm in range(numSipms):
            wfErr = np.ones(numSamples)
            #Transform both the waveform and uncertainty to mV
            wf =  wfs[ev][sipm]*range_V/np.power(2,12)*1000
            wfErr*=(range_V/np.power(2,12)*1000)/wf
            # FFT filtering rom Gary:
            denoised_waveform, freqs, fft_raw, fft_filtered = fft_denoise(wf, dt=dt, cutoff_freq=cutoff)
            wf = denoised_waveform
            wfErr*=wf
            # Do some baseline noise calculations
            baseline = np.mean(wf[:t0])
            wf = -1*(wf - baseline)
            wfEvent.append(wf)
            wfErrEvent.append(wfErr)
        wfFinal.append(wfEvent)
        wfErrFinal.append(wfErrEvent)
    return(wfFinal,wfErrFinal)

##########
################
#The intent of this function is to tag the coincident events 
#so that this could be used as a cut if desired
def tagCoincidentEvents(adjWfs):
    multSipmPulseEvents = []
    numEvents = np.shape(adjWfs)[0]
    numSipms = np.shape(adjWfs)[1]
    #numSamples = np.shape(adjWfs)[2]
    #print('tagCoincidentEvents numEvents:',numEvents)
    #print('tagCoincidentEvents numSipms:',numSipms)

    # set threshold
    peak_threshold = 5*np.std(adjWfs, axis=-1)

    for i in range(numEvents):
        tempPulseTimes = []
        tempPulseSipms = []
        for j in range(numSipms):
            peaks,_ = find_peaks(adjWfs[i][j],height=peak_threshold[i][j],distance=50,width=3)
            if(peaks.size>0 and len(tempPulseTimes)==0):
                tempPulseTimes.append(peaks)
                tempPulseSipms.append(j)
            elif(peaks.size>0 and len(tempPulseTimes)>0):
                if(np.any((tempPulseTimes >= peaks[0] - 1) & (tempPulseTimes <= peaks + 1))): 
                    tempPulseSipms.append(j)
        if(len(tempPulseSipms)>2):
            multSipmPulseEvents.append(tempPulseSipms)
        else:
            multSipmPulseEvents.append([0])
    return multSipmPulseEvents


# Put fitPulse here, but integrate prepareWaveform and droopFix
################
##########

#This is the function that does most of the work of fitting
#as well as the droop correction
def fitPulse(wf, wfErr, sample_to_us, t0):
    #First convert the time (in samples) to ns
    t_obs = np.where(wf==wf)[0]*sample_to_us
    baseline = np.mean(wf[:t0])
    baselineRMS = np.sqrt(np.mean((wf[:50])**2))
    peakBin = np.argmax(wf)
    peak = np.max(wf) #ADC value of peak
    noiseThreshold = 4.5
    initial_guess = np.array([])


    wf_corr = droopFix(wf)
   
    '''
    # Droop correction
    droop_tau = 150 # a guess
    droop_A = droop_tau * (1-np.exp(-1./droop_tau))
    S = np.zeros(np.size(wf)) # recursive correction factor
    wf_corr = np.zeros(np.size(wf))
    for i in range(np.size(wf)):
        if i < 80: # handles pre-trigger
            S[i] = 0
            wf_corr[i] = wf[i]
        elif i == 80: # acquisition t0
            S[i] = 0
            wf_corr[i] = wf[i]/droop_A
        else:
            S[i] = wf_corr[i-1] + np.exp(-1/droop_tau)*S[i-1]
            wf_corr[i] = wf[i]/droop_A + (droop_A/droop_tau)*S[i]
'''

            
    try:
        #find number of peaks. maybe need to try this, then if chi square fails, iterate about n_found
        peaks, properties = find_peaks(wf_corr,height = 6*baselineRMS, distance=15, width=10) #Ken originally had w=3
        #print('peaks, heights:',peaks, properties['peak_heights'])
        #construct initial guess and bounds arrays
        #may need to cover the scenario where there is no fit, but popt3 is held over from the previous event
        for i in range(len(peaks)):
            if i==0:
                initial_guess = [peaks[0]*sample_to_us,properties['peak_heights'][0],baseline,0.001*sample_to_us,10*sample_to_us]
                boundsLower = [peaks[0]*sample_to_us-0.25,0,-0.1,0.001*sample_to_us*0.1,10*sample_to_us*0.1]
                boundsUpper = [peaks[0]*sample_to_us+0.25,100*properties['peak_heights'][0],0.1,0.001*sample_to_us*10,10*sample_to_us*10]
            else:
                initial_guess = np.append(initial_guess, [peaks[i]*sample_to_us,properties['peak_heights'][i],baseline,0.001*sample_to_us,10*sample_to_us])
                boundsLower = np.append(boundsLower,[peaks[i]*sample_to_us-0.25,0,-0.1,0.001*sample_to_us*0.1,10*sample_to_us*0.1] )
                boundsUpper = np.append(boundsUpper,[peaks[i]*sample_to_us+0.25,100*properties['peak_heights'][i],0.1,0.001*sample_to_us*10,10*sample_to_us*10])
            bounds = np.array([boundsLower,boundsUpper])
            #print(bounds)
            # perform the fit
            popt1, pcov1 = scipy.optimize.curve_fit(newFitFunc, t_obs, wf_corr, p0=initial_guess, bounds=bounds, sigma=wfErr)
            #Let's calculate the chi-squared value for this fit
            modelY = newFitFunc(t_obs,*popt1)
            residY = np.sum(((wf_corr-modelY)/wfErr)**2)
            chiSq = residY/(len(wf)-5)

            # Manual integration - maybe we can remove this now.
            #It would be good to get a rough idea of the integral done manually
            lowerBound=peaks[0]
            upperBound=peaks[0]
            intPulse=0
            for j in range(peaks[0]-3):
                if not(wf_corr[peaks[0]-j]>wf_corr[peaks[0]-(j+1)] and wf_corr[peaks[0]-(j+1)]>wf_corr[peaks[0]-(j+2)] and wf_corr[peaks[0]-(j+2)]>wf_corr[peaks[0]-(j+3)]):
                    lowerBound=peaks[0]-(j+3)
                    #print('Found lower edge -',times[peaks[0]-(j+3)])
                    break
            for j in range(len(wf_corr)-(peaks[0]+4)):
                if not(wf_corr[peaks[0]+j]>wf_corr[peaks[0]+(j+1)] and wf_corr[peaks[0]+(j+1)]>wf_corr[peaks[0]+(j+2)] and wf_corr[peaks[0]+(j+2)]>wf_corr[peaks[0]+(j+3)]):
                    upperBound=peaks[0]+(j+3)
                    #print('Found upper edge -',times[peaks[0]+(j+3)])
                    break
            for j in range(upperBound-lowerBound):
                #print(wf_corr[lowerBound+j])
                intPulse+=wf_corr[lowerBound+j]
            return popt1, pcov1, wf_corr, wfErr, t_obs, chiSq, intPulse*sample_to_us
        else:
            return [-1000,-1000,-1000,-1000,-1000], [np.inf,np.inf,np.inf],wf_corr, t_obs, -100, 100, -100
    except RuntimeError:
            return [-1000,-1000,-1000,-1000,-1000], [np.inf,np.inf,np.inf],wf_corr, t_obs, -100, 100, -100
    except ZeroDivisionError:
            return [-1000,-1000,-1000,-1000,-1000], [np.inf,np.inf,np.inf],wf_corr, t_obs, -100, 100, -100

def getFitValues(ev, numEvToDo=-1, t0=80):
    if ev is None:
        return {
          "thit": np.array([]),
          "area": np.array([]),
          "baseline": np.array([]),
          "fall_time": np.array([]),
          "rise_time": np.array([]),
          "chi2": np.array([]),
          "evnum": np.array([]),
          "coinc": np.array([])
        }


    wfs = ev["scintillation"]["Waveforms"]

    try: # try loading decimation, otherwise assume it is 1
        decimation = ev["run_control"]["caen"]["global"]["decimation"]
    except:
        decimation = 1
    sample_rate =  SAMPLE_FREQ/(2**decimation)
    sample_to_us = 1/sample_rate
    
    #The array is [#events,#sipms,#samples]
    numEvents = np.shape(wfs)[0]
    numSipms = np.shape(wfs)[1]
    numSamples = np.shape(wfs)[2]

    if numEvToDo > 0 and numEvToDo < numEvents:
        numEvents = numEvToDo

    #print('getFitValues numEvents:',numEvents)
    #print('getFitValues numSipms:',numSipms)
    #print('getFitValues numsSamples:',numSamples)

    #First let's remove the baseline and flip the pulses
    adjWfs,adjWfsErr = prepareWaveforms(wfs, t0)

    #Next let's do the simple test to sort out coincident events
    coincEvents = tagCoincidentEvents(adjWfs)

    #Then let's get down to finding the fit values
    #Start by setting up the values that we want to store
    sipmsFitT0 = np.full((numSipms, numEvents), np.nan)
    sipmsFitArea = np.full((numSipms, numEvents), np.nan)
    sipmsFitBaseline = np.full((numSipms, numEvents), np.nan)
    sipmsFitFallTime = np.full((numSipms, numEvents), np.nan)
    sipmsFitRiseTime = np.full((numSipms, numEvents), np.nan)
    sipmsManual = np.full((numSipms, numEvents), np.nan)
    sipmsChiSq = np.full((numSipms, numEvents), np.nan)
    sipmsEvNum = np.full((numSipms, numEvents), np.nan)
    sipmsCoinc = np.full((numSipms, numEvents), np.nan)

    #Now that those are set up, let's actually do this thing
    for i in range(numEvents):
        for j in range(numSipms):
            fitParams, fitCov, denoisedWf, wfErr, times, chiSq, intPulse = fitPulse(adjWfs[i][j], adjWfsErr[i][j], sample_to_us, t0)
            sipmsFitT0[j][i] = fitParams[0]
            sipmsFitArea[j][i] = fitParams[1]
            sipmsFitBaseline[j][i] = fitParams[2]
            sipmsFitFallTime[j][i] = fitParams[3]
            sipmsFitRiseTime[j][i] = fitParams[4]
            sipmsManual[j][i] = intPulse
            sipmsChiSq[j][i] = chiSq
            sipmsEvNum[j][i] = i
            if j in coincEvents[i]:
                sipmsCoinc[j][i] = 1
            else:
                sipmsCoinc[j][i] = 0

    return {
      "thit": np.array(sipmsFitT0),
      "area": np.array(sipmsFitArea),
      "baseline": np.array(sipmsFitBaseline),
      "fall_time": np.array(sipmsFitFallTime),
      "rise_time": np.array(sipmsFitRiseTime),
      "chi2": np.array(sipmsChiSq),
      "evnum": np.array(sipmsEvNum),
      "coinc": np.array(sipmsCoinc)
    }

def SiPMFitterBatched(ev, t0=80, nwvf_batch=1000, maxwvf=0, progress=False, njob=1):
    return BatchSiPMs.BatchSiPMs(ev, getFitValues, t0=t0, 
        nwvf_batch=nwvf_batch, maxwvf=maxwvf, progress=progress, njob=njob)


if __name__ == "__main__":
    from GetEvent import GetEvent

    TEST_RUN = "/exp/e961/data/SBC-25-daqdata/20250703_11/"
    TEST_EVENT = 0
    out = getFitValues(GetEvent(TEST_RUN, TEST_EVENT, strictMode=False, lazy_load_scintillation=False), 1)
    print(out["thit"])

