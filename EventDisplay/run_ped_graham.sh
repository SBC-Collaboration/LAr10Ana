#!/bin/bash

# Note: as this requires a graphical session, the recommended method is to connect with 
#   TigerVNC Viewer to gra-vdi.computecanada.ca with the ~/.vnc/x509_ca.pem SSL certificate 
#   symlinked as directed by https://docs.computecanada.ca/wiki/VNC (if necessary for your OS).
#   This script should then be run from a terminal in that VNC session.
#
#   For convenience, you may symlink this to your home directory, e.g.
#     $ cd; ln -s projects/rrg-kenclark/pico/PEDsvn/run_ped_graham.py


# Load the required python-related modules for use on CC/graham's graphical VDI login nodes
#  (note: one could run all of this interactively, but modules unload when a shell terminates)
# module load nixpkgs/16.09 python/3.7.0 scipy-stack/2019a
echo "Loading python and scipy-stack"
module load python/3.7.0 scipy-stack/2019a

# Then activate a virtualenv in which we've `pip install`-ed Pillow (see: migration_notes.txt)
echo "Activating venv"
# source /project/6007972/pico/venv_ped_2/bin/activate
source /project/6007972/pico/venv_ped_3/bin/activate

# Now we can actually run PED itself from a graphical session
echo "Attempting to run Event Display"
python /project/6007972/pico/EventDisplay/eventdisplay/ped.py

# Deactivate the virtualenv
deactivate
