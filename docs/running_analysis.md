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
This is the script to run to submit multiple jobs. It checks the folder of both the raw data tarball and the analysis outputs. For each run, it checks if the recon output doesn't exist, or current repository tag (like `v0.0.7`) is newer than the existing output for this run. If either is true, it will proceed to submitting a job for this run. It keeps track of the run ID and job ID of submitted jobs at `~/.cache/sbc/job_list.csv`, and also uses a lockfile to prevent multiple instances running at the same time.  

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

### 4. `clean_up.sh`
This is another piece of script that needs to be run explicitly. Since all data and output files are in the scratch part of PNFS, they are not persistent. For each run in the `grid_output` folder, it checks the exit code of the log file. If the code is not zero, then the output folder is deleted. If it's zero, meaning success, it will then compare the version of the output to the run existing in the destination (`/exp/e961/data/SBC-25-recon/dev-output` in most cases). If the version is the same or newer, then the output in the destination will be overwritten by the new output. If not, the new output will be discarded. Then, it removes the run from `~/.cache/sbc/job_list.csv`. The raw data tarball is left in temp_data folder, to avoid repeated file operations. 
It has a tag `--tag` or `-t`. If a version (`v0.1.0` for example) is passed with the tag, then it will look instead at the `grid_output_v0.1.0` folder, and move all finished data back to `/exp/e961/data/SBC-25-recon/v0.1.0`.

### 5. `batch_clean_up.sh`
This is basically a wrapper script that detects all the versions available in the `/pnfs/coupp/scratch/coupppro` folder, and runs `clean_up.sh` on each of those versions, including the main one without a version. Then it cleans up any empty folders in the temp folder.

### Cron jobs
The analysis is run on any new data, and on all data when a new version is added to `~/LAr10Ana` in the `coupppro` user home directory. This is done by cron jobs. Every six hours, it runs `batch_run_gridjob.sh` to submit new jobs. Every hour, it runs `batch_clean_up.sh` to move any completed analysis to the destination. The log file of each of the scrips are saved to `~/.cache/sbc/`. To see the current cron job of a user, run `crontab -l`. To edit, run `crontab -e`.

## Deploying a new version
Here's the recommended steps when running a new version on the grid.
1. Merge all changes that should be included in the new version into the main branch. 
2. Test the code works, and can produce the expected recon output files. On the GPVM or EAF, pull the latest version of main branch, then edit the last section of the `grid_jobs/EventDealer.py` to a run that has all data streams to be tested, the output path to a personal temp folder, and the updated module list. Then run
    ```
    source setup.py  # set up conda environment
    python grid_jobs/EventDealer.py
    ```
    and check that the analysis modules finish for this run, and the output files are as expected.
3. Once you are happy with the code, bump the version number in `docs/conf.py` (for the Sphinx documentation), and update change log at `docs/changelog.md`. Push the commit.
4. Tag this commit with the new version number. To do it locally, run
    ```
    git tag "vx.y.z"
    git push --tags
    ```
    or do it on GitHub by creating a new release. 
5. Update the `LAr10Ana` repository on coupppro@couppsbcgpvm01.fnal.gov to pull the new tag and commit. 
6. It is ready to be deployed. Either wait for the cron job to run, or manually start by running 
    ```
    grid_jobs/batch_run_grid_jobs.sh
    ```
    To run the specific version, run 
    ```
    grid_jobs/batch_run_grid_jobs.sh -t vx.y.z
    ```

## Adding a new module
To add a new module to run by EventDealer on the grid, following these steps.
1. Create a feature branch based on the latest main branch. Ideally give it a descriptive name.
2. If you are working on this branch for more than a couple days, consider:
    - Commit any work along the way. Push your changes to github, the new branch should show up.
    - Rebase onto the main branch periodically when new commits are available in the main branch. This means syncing the new commits from the main branch before the first unique branch in the feature branch.
3. Add the new module as a separate python file inthe `ana/` directory. See existing modules for guideline of structure. Mainly, each module should have a main function that takes the output of `GetEvent` as an argument, and returns a dictionary of pre-defined keys and values. Every value should be a list / array in the dictionary should have the same size, but that size can be different for each event. This dictionary does not need to include the run and event ID, as that will be added by the Event Dealer.
4. Edit `ana/EventDealer.py` to add the following lines. For example, this uses the `BubbleFinder` module.
    - At the top, import the main function from the module
    ```
    from ana.BubbleFinder import BubbleFinder as bf
    ```
    - Include the module in the list os analyses to run
    ```
        ANALYSES = {
        ...
        "bubble": bf,
        ...
    }
    ```
    - Add to the list of checks the condition to skip this analysis module. I.e. if an event does not have camera images, don't run the bubble finder module
    ```
    ...
    elif p == "bubble" and not data["cam"]["loaded"]:
    print(f"Skipping {p} analysis -- camera data not loaded.")
    continue
    ...
    ```
    - Add the module into the list of analyses to run in the main function
    ```
    if __name__ == "__main__":
    if len(sys.argv) > 1:
        ProcessSingleRun(
            rundir=sys.argv[1],
            recondir=sys.argv[2],
            process_list = [..., "bubble", ...])
    else:
        ProcessSingleRun(
            rundir="/exp/e961/data/SBC-25-daqdata/20260212_1.tar",
            recondir="/PATH/TO/TEST/OUTPUT/DIRECTORY/"
            process_list = [..., "bubble", ...])
    ```
    - Then test run the new version of Event Dealer on GPVM, to make sure it can handle the output format of the new module without throwing any errors.
    - Add documentation of the data output to `docs/data_format.md`. 
    - Commit and push changes to the feature branch, open a pull request to main branch.
    - When ready, create a new version (instruction above), and start running on the grid~