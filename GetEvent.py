import os
import numpy as np
from sbcbinaryformat import Streamer, Writer
from PIL import Image
import json
import warnings

full_loadlist = [
    "acoustics",
    "scintillation",
    "cam",
    "event_info",
    "plc",
    "slow_daq",
    "run_info",
    "run_control"
]

def GetScint(ev, start=None, end=None, length=None):
    out_ev = dict([(k, v.copy()) for (k, v) in ev.items()]) # copy input

    for key in ev["scintillation"].keys():
        if key == "loaded" or key == "length" or key == "sample_rate" or key == "EventCounter": # skip helper keys
            continue

        out_ev["scintillation"][key] = ev["scintillation"][key](start=start, end=end, length=length)

    return out_ev

def GetEvent(rundirectory, ev, *loadlist, strictMode=True, lazy_load_scintillation=True):
    event = dict()
    event_dir = os.path.join(rundirectory, str(ev)) 

    for key in full_loadlist:
        event[key] = dict(loaded=False)

    if len(loadlist) == 0:
        loadlist = full_loadlist
    elif loadlist[0][0] == "~":
        loadlist = [l for l in full_loadlist if l not in [s.lstrip("~") for s in loadlist]]

    if "acoustics" in loadlist:
        acoustic_file = None
        for fname in os.listdir(event_dir):
            if fname.startswith("acoustics"):
                acoustic_file = os.path.join(event_dir, fname)
                break

        if acoustic_file is None:
            if strictMode: 
                raise FileNotFoundError("No acoustics file present in the run directory. To disable this error, either pass strictMode=False, or remove 'acoustics' from the loadlist")
            else:
                warnings.warn("No acoustics file present in the run directory. Data will not be available in the returned dictionary.")
        else:
            acoustic_data = Streamer(acoustic_file).to_dict()
            event["acoustics"]["loaded"] = True
            for k, v in acoustic_data.items():
                event["acoustics"][k] = v
        
    if "scintillation" in loadlist:
        scint_file = os.path.join(event_dir, "scintillation.sbc")

        if not os.path.exists(scint_file):
            if strictMode: 
                raise FileNotFoundError("No scintillation file present in the run directory. To disable this error, either pass strictMode=False, or remove 'scintillation' from the loadlist")
            else:
                warnings.warn("No scintillation file present in the run directory. Data will not be available in the returned dictionary.")
        else:
            scint = Streamer(scint_file)
            if lazy_load_scintillation:
                for c in scint.columns:
                    event["scintillation"][c] = lambda start=None, end=None, length=None: scint.to_dict(start=start, end=end, length=length)[c]
                event["scintillation"]["length"] = scint.num_elems
            else:
                scint = scint_data.items()
                for k, v in scint_data.items():
                    event["scintillation"][k] = v

            event["scintillation"]["loaded"] = True

    if "cam" in loadlist:
        event["cam"]["loaded"] = True
        for cam_ind in range(1, 4):
            event["cam"]["c%i" % cam_ind] = {}
            cam_file = os.path.join(event_dir, "cam%i-info.csv" % cam_ind)
            if not os.path.exists(cam_file):
                if strictMode: 
                    raise FileNotFoundError("Missing camera file in the run directory. To disable this error, either pass strictMode=False, or remove 'cam' from the loadlist")
                else:
                    warnings.warn("Missing camera file in the run directory. Data will not be available in the returned dictionary.")
                continue

            cam_data = np.loadtxt(cam_file, delimiter=",", skiprows=1)
            cam_data_headers = ["index"]
            with open(cam_file) as f:
                first_line = f.readline()
                cam_data_headers += [s for s in first_line.rstrip("\n").split(",") if s]
            for h, d in zip(cam_data_headers, cam_data):
                event["cam"]["c%i" % cam_ind][h] = d

        for fname in os.listdir(event_dir):
            if fname.startswith("cam") and fname.endswith(".png"):
                img_file = os.path.join(event_dir, fname)
                cam_ind = int(fname[3])
                frame_ind = int(fname[8:10])
                event["cam"]["c%i" % cam_ind]["frame%i" % frame_ind] = np.array(Image.open(img_file).convert("RGB"))

    if "event_info" in loadlist:
        event_file = os.path.join(event_dir, "event_info.sbc")

        if not os.path.exists(event_file):
            if strictMode: 
                raise FileNotFoundError("No event_info file present in the run directory. To disable this error, either pass strictMode=False, or remove 'event_info' from the loadlist")
            else:
                warnings.warn("No event_info file present in the run directory. Data will not be available in the returned dictionary.")
        else:
            event_data = Streamer(event_file).to_dict()
            event["event_info"]["loaded"] = True
            for k, v in event_data.items():
                event["event_info"][k] = v

    if "slow_daq" in loadlist:
        slow_daq_file = os.path.join(event_dir, "slow_daq.sbc")
        if not os.path.exists(slow_daq_file):
            if strictMode: 
                raise FileNotFoundError("No slow_daq file present in the run directory. To disable this error, either pass strictMode=False, or remove 'slow_daq' from the loadlist")
            else:
                warnings.warn("No slow_daq file present in the run directory. Data will not be available in the returned dictionary.")
        else:
            slow_daq_data = Streamer(slow_daq_file).to_dict()
            event["slow_daq"]["loaded"] = True
            for k, v in slow_daq_data.items():
                 event["slow_daq"][k] = v

    if "plc" in loadlist:
        plc_file = os.path.join(event_dir, "plc.sbc")
        if not os.path.exists(plc_file):
            if strictMode: 
                raise FileNotFoundError("No plc file present in the run directory. To disable this error, either pass strictMode=False, or remove 'plc' from the loadlist")
            else:
                warnings.warn("No plc file present in the run directory. Data will not be available in the returned dictionary.")
        else:
            plc_data = Streamer(plc_file).to_dict()
            event["plc"]["loaded"] = True
            for k, v in plc_data.items():
                event["plc"][k] = v

    if "run_info" in loadlist:
        run_info_file = os.path.join(rundirectory, "run_info.sbc")
        if not os.path.exists(run_info_file):
            if strictMode: 
                raise FileNotFoundError("No run_info file present in the run directory. To disable this error, either pass strictMode=False, or remove 'run_info' from the loadlist")
            else:
                warnings.warn("No run_info file present in the run directory. Data will not be available in the returned dictionary.")
        else:
            run_info_data = Streamer(run_info_file).to_dict()
            event["run_info"]["loaded"] = True
            for k, v in run_info_data.items():
                event["run_info"][k] = v

    if "run_control" in loadlist:
        run_ctrl_file = os.path.join(rundirectory, "rc.json")
        if not os.path.exists(run_ctrl_file):
            if strictMode: 
                raise FileNotFoundError("No run_control file present in the run directory. To disable this error, either pass strictMode=False, or remove 'run_control' from the loadlist")
            else:
                warnings.warn("No run_control file present in the run directory. Data will not be available in the returned dictionary.")
        else:
            with open(run_ctrl_file, "r") as f:
                run_ctrl_data = json.load(f)
            event["run_control"]["loaded"] = True
            for k, v in run_ctrl_data.items():
                event["run_control"][k] = v
            sample_rate_str = event['run_control']['acous']['sample_rate'].strip().upper()

            if "MS/S" in sample_rate_str:
                sample_rate = int(sample_rate_str.replace("MS/S", "").strip()) * 1_000_000
            elif "KS/S" in sample_rate_str:
                sample_rate = int(sample_rate_str.replace("KS/S", "").strip()) * 1_000
            elif "S/S" in sample_rate_str:
                sample_rate = int(sample_rate_str.replace("S/S", "").strip())
            else:
                raise ValueError(f"Unrecognized sample rate format: '{sample_rate_str}'")

            decimation = event['run_control']['caen']['global']['decimation']

            event["acoustics"]["sample_rate"] = sample_rate

            if "scintillation" in loadlist:
                event['scintillation']['sample_rate'] = 62500000 / (2**decimation)

    return event
