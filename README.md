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

# Running Analysis on the FermiGrid

TODO
