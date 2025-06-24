# Getting started on FNAL servers

1. Send a ticket requesting access to the coupp organization with the [Fermilab Service Desk](https://fermi.servicenowservices.com/wp)
2. Once you have access, ssh into the general purpose virtual machine (GPVM)
   ```
   kinit <username>@FNAL.GOV
   ssh -KYX <username>@couppsbcgpvm01.fnal.gov
   ```
3. Navigate to the app directory area, then make a new directory for yourself
   ```
   cd /exp/e961/app/users
   mkdir <username>
   cd <username>
   ```
4. Get the code! Setup an [SSH key for github](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
   if you do not already have one (follow the "Generating a new SSH key" steps), and then [add it to you github account](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account).
   Also copy the key to your home directory on `couppsbcgpvm01.fnal.gov` to: `~/.ssh/`.
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
kinit <username>@FNAL.GOV
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

The first time you run a notebook, you'll want to add the conda environment to your list of available jupyterlab kernels. To do this, navigate to your `LAr10Ana` directory, and run the commands:
```
cd /path/to/LAr10Ana
source setup.sh
source jupyter_init.sh
```

Now you're ready to start a jupyter notebook! To do this, start from a terminal inside `LAr10Ana` where you have already activated the conda environment (by running `setup.sh`). Then, run:
```
jupyter notebook --no-browser --port <port>
```

Where `<port>` is the port you selected in the ssh. A link will print in the terminal. Open it in the browser to access the jupyter notebook. For some example notebooks for 
accessing SBC-LAr10 data, navigate to `nb/examples`. Make sure to select the `SBC conda (env)` kernel when starting a new notebook.

# Running Analysis on the FermiGrid

TODO
