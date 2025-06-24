import numpy as np

import warnings
# np.nanmin will raise a RuntimeWarning if all inputs are nan (i.e., no hit in a readout). Ignore these
warnings.filterwarnings("ignore", category=RuntimeWarning)

def PhotonT0(ev):
    default_out = dict(
        t0=np.array([]),
        amp=np.array([]),
    )
    out = default_out

    # SiPMPulses should've been run previously on the event.
    # This module combines the hit data from that module to get a final T0 and voltage

    hit_times = ev["hit_t0"]
    hit_ampls = ev["hit_amp"]

    # earliest hit time is the t0
    # Ignore nans, which corresponds to no hit
    combined_t0 = np.nanmin(hit_times, axis=0)

    # Sum of amplitudes is the total amplitude
    # Ignore nans, which corresponds to no hit
    combined_amp = np.nansum(hit_ampls, axis=0)

    out["t0"] = combined_t0
    out["amp"] = combined_amp

    return out

if __name__ == "__main__":
    import os
    from DataHandling import ReadBinary, GetSBCEvent
    from SiPMPulses import SiPMPulses

    path = "/bluearc/storage/SBC-22-caendata"
    file = "202206201045.bin"
    
    filename = os.path.join(path, file)
    
    data = ReadBinary.ReadBlock(filename)
    sipm_data = SiPMPulses(data)
    print(PhotonT0(sipm_data))
