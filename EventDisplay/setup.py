#!/usr/bin/env python

from setuptools import setup

with open('requirements.txt') as f:
    requirements = f.readlines()

setup(
    name='EventDisplay',
    description='Display for PICO events',
    author='PICO collaboration',
    packages=[
	'eventdisplay'
    ],
#    dependency_links = [
#        'https://github.com/picoexperiment/PICOcode/tarball/make-library#egg=PICOcode=0.0.1'
#    ],
    install_requires=requirements
)
