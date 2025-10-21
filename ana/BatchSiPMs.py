# from tqdm.auto import tqdm as tq
from tqdm import tqdm 
from multiprocess import Pool
from GetEvent import GetScint
import numpy as np

def BatchSiPMs(ev, ana_f, nwvf_batch=1000, maxwvf=0, progress=False, njob=1, **f_kwargs):
    # load defaults
    output = ana_f(None)

    if ev is None:
        return output

    # Call SiPMPulses in batches
    nwvf = ev["scintillation"]["length"]
    if nwvf > maxwvf:
        nwvf = maxwvf
    print("BATCHING %i pulses" % nwvf)

    batched_outputs = []

    itr = list(range(0, nwvf, nwvf_batch))
    nitr = len(itr)

    def applyf(start):
        end = min(start + nwvf_batch, nwvf)
        return ana_f(GetScint(ev, start=start, end=end), **f_kwargs)

    if njob > 1:
        pool = Pool(processes=njob)
        itr = pool.imap_unordered(applyf, itr)
    else:
        itr = map(applyf, itr)

    if progress:
        itr = tqdm(itr, total=nitr)

    for pulses in itr:
        batched_outputs.append(pulses)

    if njob > 1:
        pool.close()

    # concatenate the outputs
    for key in output.keys():
        output[key] = np.concatenate([b[key] for b in batched_outputs], axis=1) 
        
    return output
