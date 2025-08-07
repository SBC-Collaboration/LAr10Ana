This is a miniconda, SBC event display, or PED, designed for speedy event-by-event navigation through raw data (camera images, piezo traces, temperature/pressure data) with optional ability of doing hand scanning. When reconstructed data is available, PED, also does speedy cut-based navigation and display of reconstructed data (bubble position, acoustic parameter).

Depdencies are managed in the parent directory in `conda_init.sh` and `setup.sh`

**Notes:**
The config files allow specification of raw/reconstructed data location for a given dataset, as well as default image orientation, number of cameras, trigger frame number and recorded frame range, default piezo and PLC temperature variable names.

**IMPORTANT:**
In order for PED to read raw and/or reconstructed data (i.e. in order for it to work), a "convert" script needs to be run to create an .npy file (one for raw, one optionally for reco) for each dataset of interest. For raw data this traverses the folder-based hierarchical "run/eventNum" data structure and compactly stores the navigable information for fast retrieval. For reco data the script reads the merged_all.txt file and stores it in a speedier format for fast retrieval. There are a number of convert scripts. The most recent, designed for ComputeCanada operation, are the following, to be run in this order:

- convert_raw_to_npy_run_by_run.py (this creates npy files for each raw run that doesn't already have one, untarring if necessary)
- merge_raw_run_npy.py (this merges the npy files for individual runs into a single raw_events.npy file)
- convert_reco_to_npy_and_reindex_raw_npy.py (this creates reco_events_all.npy from merged_all.txt and re-indexes the raw_events.npy file for making fast cuts on reco data, and then makes a culled reco_events.npy file only containing the reco data for the raw data that is present)

The `update_npy_data.sh` script is currently running as a cronjob on the coupp
server. This can also be run manually on coupp to rebuild `raw_events.npy`
