#!/bin/bash

# Commands to create a new conda environment with the packages needed for SBC LAr10 analysis

# NOTE: YOU PROBABLY DON'T NEED TO RUN THIS. USER CODE SHOULD BE RUN WITH THE CENTRAL CONDA ENVIRONMENT.

conda create --name env python jupyterlab notebook numpy pandas scipy matplotlib pip nbstripout opencv pymysql spyder plotly ipympl tqdm multiprocess pywavelets scikit-image diplib
conda activate env
pip install git+https://github.com/SBC-Collaboration/SBCBinaryFormat.git@master#subdirectory=python
