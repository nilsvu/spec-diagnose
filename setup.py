#!/usr/bin/env python

from distutils.core import setup

setup(
    name='spec_diagnose',
    version='0.0.1',
    description="Utility for diagnosing SpEC simulations",
    author="SXS collaboration",
    url="https://black-holes.org",
    packages=['spec_diagnose'],
    install_requires=['h5py', 'matplotlib', 'numpy', 'tqdm'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
