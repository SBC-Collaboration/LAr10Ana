# Getting started on FNAL servers

1. Send a ticket requesting access to the COUUP affiliation with the [Fermilab Service Desk](https://fermi.servicenowservices.com/com.glideapp.servicecatalog_cat_item_view.do?v=1&syspa[…]alog_view=catalog_default&sysparm_view=catalog_default). This can take up to 2 business days.
2. Once you have access, ssh into the general purpose virtual machine (GPVM)
   ```
   kinit <username>@FNAL.GOV
   ssh -KYX <username>@couppsbcgpvm01.fnal.gov
   ```
   If you see errors relating to permissions in your home directory, try using 
   `kinit -f <username>>@fnal.gov` to make sure your credentials are being forwarded.
3. Navigate to the app directory area, then make a new directory for yourself
   ```
   cd /exp/e961/app/users
   mkdir <username>
   cd <username>
   ```
4. Get the code! Setup an [SSH key for github](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
   if you do not already have one (follow the "Generating a new SSH key" steps), and then [add it to you github account](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).
   Also copy the key to your home directory on `couppsbcgpvm01.fnal.gov` to: `~/.ssh/`. If you ran the ssh keygen on the `couppsbcgpvm01.fnal.gov` machine, it will be there already.
   Then clone the directory into you app users area:
   ```
   git clone git@github.com:SBC-Collaboration/LAr10Ana.git
   ```
   Now you have access to the code on the Fermilab servers! See below for instructions on other things you may want to do next:
   how to setup the SBCLAr10 conda environment, how to log back into the fermilab servers, how to submit grid jobs, and others.

# Logging back into the Fermilab servers

Once you have your account setup, to get back to your app area from your laptop, run:

```
ssh -KYX <username>@couppsbcgpvm01.fnal.gov
cd /exp/e961/app/users/<username>
```

You may get an authentication error if you have not run `kinit` in about a day. In that case, before `ssh`-ing, run:

```
kinit -f <username>@FNAL.GOV
```

# Activating the SBCLAr10 Conda Environment

SBC analysis mostly relies on python code, with a few external dependencies. The code is managed in a central conda environment.
Once you have the LAr10Ana repository downloaded, navigate to `LAr10Ana/`. Then run:

```
source setup.sh
```

This will setup the conda environment in your current terminal.

# Running jupyter notebooks

To run a jupyter notebook, you'll want to start a server from your terminal and then open it in a browser. First, login to the Fermilab server and forward port `8888`:

```
ssh -KYX <username>@couppsbcgpvm01.fnal.gov -L 8888:localhost:8888
```

If you see an error in the login that port 8888 is taken, then try again with a different port (8889, for example).

The first time you run a notebook, you'll want to add the conda environment to your list of available jupyterlab kernels. To do this, navigate to your `LAr10Ana` directory (`cd /path/to/LAr10Ana`), and run the commands:

```
source setup.sh # This activates the conda environment. You'll need to run it every time you create a new terminal connection.
source jupyter_init.sh # This tells jupyter about your conda environment. You only need to run it the first time you run the notebook.
```

Now you're ready to start a jupyter notebook! To do this, start from a terminal inside `LAr10Ana` where you have already activated the conda environment (by running `setup.sh`). Then, run:

```
jupyter notebook --no-browser --port <port>
```

Where `<port>` is the port you selected in the ssh. A link will print in the terminal. Open it in the browser to access the jupyter notebook. For some example notebooks for
accessing SBC-LAr10 data, navigate to `nb/examples`. Make sure to select the `SBC conda (env)` kernel when starting a new notebook.

# Running Event Display

To run the event display, you'll need to launch a VNC session on a free port.

## Check for VNC ports

See which ports are already in use:

```
/exp/e961/app/home/coupp/bin/RegisterVNCDesktop
```

Pick a two-digit number that isn't taken. For example, 50 would mean your VNC server runs on port 5950.

Register your choice:

```
/exp/e961/app/home/coupp/bin/RegisterVNCDesktop 50 "Optional short description"
```

## Setup `xstartup`

If you don’t already have ~/.vnc/xstartup, make one. Also create ~/.vnc if needed.

Here's a working example:

```#!/bin/sh

# Unset conflicting environment variables
unset SESSION_MANAGER
unset DBUS_SESSION_BUS_ADDRESS
unset PYTHONHOME PYTHONPATH LD_LIBRARY_PATH LD_PRELOAD

# Load X resources
[ -f $HOME/.Xresources ] && xrdb $HOME/.Xresources

# Start IceWM
exec icewm
```

Make sure it’s executable:

```
chmod +x ~/.vnc/xstartup
```

## Start the VNC Server

Use the display number you picked earlier:

```
vncserver :PORT_NUMBER
```

## Run the Event Display

You will need to port forward your local VNC port to your localhost machine. This
can be done using `ssh`

`ssh -L <local_port>:localhost:<remote_port> <user>@<remote_host>`

After tunneling you can run the following commands to run the event display

```
export DISPLAY=:PORT_NUMBER
cd /exp/e961/app/LAr10Ana
source setup.sh
python EventDisplay/eventdisplay/ped.py
```

## VNC Software

You will need to download VNC software and connect to `localhost:PORT_NUMBER` to
view the event display.

- For MacOS we recommend using the built in "Screen Sharing" software.
- For Windows we recommend using TigerVNC
- For Chromebook, try using Remmina (Linux app)

# Running Analysis 
`EventDealer` is the piece of code that runs a list of analysis modules in the `LAr10Ana/ana/` folder on a single run, and saves the all output to file. You can run it directly to test and debug, or to run it on the FermiGrid using the scripts on all data.

## Running Directly (on GPVM or EAF)
To test if an analysis module is working properly, you can run the `EventDealer` on a single run on `couppsbcgpvm01.fnal.gov` or [EAF](https://analytics-hub.fnal.gov). You'll get a more complete error message in the output, and don't need to wait for the job to be scheduled. The file IO path will be different from when running on the grid, so it doesn't help with debugging permission issues, etc.
1. Clone the LAr10Ana repository in your home directory or a personal working directory like `/exp/e961/app/users/<username>`.
2. Edit `LAr10Ana/grid_jobs/EventDealer.py`, change the path to the run data tarball, the path to save the output, and the processes to run.
```
if __name__ == "__main__":
    if len(sys.argv) > 1:
        ProcessSingleRun(
            rundir=sys.argv[1],
            recondir=sys.argv[2],
            process_list = ["run", "event", "exposure", "scintillation", "scint_rate", "bubble", "reco", "pressure_t0", "t_expansion"])
    else:
        ProcessSingleRun(
            rundir="/exp/e961/data/SBC-25-daqdata/20260221_0.tar",
            recondir="/nashome/z/zsheng/test", # Use your own directory for testing~
            process_list = ["run", "event", "exposure", "scintillation", "scint_rate", "bubble", "reco", "pressure_t0", "t_expansion"])
```
3. Activate conda environment
```
source LAr10Ana/setup.sh
```
4. Run EventDealer in one of two ways
```
# using info in the else clause
python LAr10Ana/grid_jobs/EventDealer.py

# provide info in argument
python LAr10Ana/grid_jobs/EventDealer.py /exp/e961/data/SBC-25-daqdata/20260221_0.tar /nashome/z/zsheng/test
```

## Running on the FermiGrid
The `LAr10Ana/grid_jobs/` directory also contains a few bash scripts to help with submitting and running jobs on the Fermi grid. The jobs need to be submitted from `couppsbcgpvm01.fnal.gov`. For production mode and automated submission, use `coupppro` user, which is set up with managed kerberos ticket.

### 1. `batch_run_gridjob.sh`
This is the script to run to submit multiple jobs. It checks the folder of both the raw data tarball and the analysis outputs. For each run, it checks if the recon output doesn't exist, or current repository tag (like `v0.0.7`) is newer than the existing output for this run. If either is true, it will proceed to submitting a job for this run. It keeps track of the run ID and job ID of submitted jobs at `~/.cache/sbc_job_list.csv`, and also uses a lockfile to prevent multiple instances running at the same time.  

Three tags are available for this script.  
- `--verbose` or `-v`: Verbose mode, more detailed logging. Passed to `run_gridjob.sh`.
- `--force-rerun` or `--force` or `-f`: Force mode, submit jobs for all runs without the checks
- `--production` or `--pro` or `-p`: Production mode, for automated submission using production role permission and a set of different paths. Passed to `run_gridjob.sh`.

### 2. `run_gridjob.sh`
This script is used to submit jobs for a single run. It can be used indepently or by `batch_run_gridjob.sh`. It will copy the raw data tarball to the PNFS where the grid nodes have access to, make a tarball of the LAr10Ana repository, create output path in PNFS if it doesn't yet exist, and actually submit the job.  

The working directory for each user is at `/pnfs/coupp/scratch/users/${USER}` to allow proper permission. Two folders will be created inside: `temp_data` to store a copy of the data tarball, and `grid_output` to store the output files of the analysis.  

When submitting jobs using the [jobsub_lite](https://github.com/fermitools/jobsub_lite) program, it also calcuates the appropriate disk space, RAM, and run time to request. Too much resources requested may mean longer wait time before the job is scheduled, while too little requested may result in job exceeding limit and being held. The resoruces are calculated based on the size of the tarball. This can be improved to distinguish between runs with many events, vs runs with few events but large files in each (like long scintillation runs).
- Disk: 3x size + 5GB, min 5GB, max 500GB
- RAM: 2x size + 2GB, minimum 2GB, max 16GB
- Run time: 1h/5GB x size + 2h, minimum 2h, maximum 24h  

Two tags are available:
- `--verbose` or `-v`: Verbose mode, more detailed logging. 
- `--production` or `--pro` or `-p`: Production mode. This only works if submitted as `coupppro@couppsbcgpvm01.fnal.gov`. The working directory in production mode is at `/pnfs/coupp/scratch/coupppro` instead. The grid jobs are submitted with production role using the managed kerberos ticket. (This presumply allows for priority and better resources allocation compared to normal analysis role, but I haven't really noticed a difference, since not a lot of people in the coupp group are running jobs. -ZS)  

### 3. `gridjob.sh`
This script is the one that actually runs on the grid node. It configures conda environment, copies the data from PNFS to the internal storage of the node, untars it, and runs `EventDealer.py`. After `EventDealer.py` is finished, is copies the output folder from the node back to the PNFS.

### 4. `cleanup.sh`
This is another piece of script that needs to be run explicitly. Since all data and output files are in the scratch part of PNFS, they are not persistent. For each run in the `grid_output` folder, it checks the exit code of the log file. If the code is not zero, then the output folder is deleted. If it's zero, meaning success, it will then compare the version of the output to the run existing in the destination (`/exp/e961/data/SBC-25-recon/dev-outputs` in most cases). If the version is the same or newer, then the output in the destination will be overwritten by the new output. If not, the new output will be discarded. Then, it removes the run from `~/.cache/sbc_job_list.csv`. The raw data tarball is left in temp_data folder, to avoid repeated file operations. 

## Further Resources for Documentation

The [Run Control documentation](https://runcontrol.readthedocs.io/latest/) has information on the SBC Run Control, DAQ system, and raw data varaiables.

A sub page in this documentation has information on the reconstructed/analysis output file variables: https://runcontrol.readthedocs.io/latest/recon_data_format.html.
