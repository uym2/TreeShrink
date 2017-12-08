#!/usr/bin/env python

#############################################################################
##  this file is part of TreeShrink.
##  see "license.txt" for terms and conditions of usage.
#############################################################################


"""
Package setup and installation.
"""

import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages

from datetime import datetime
import os
import platform
import sys
import treeshrink

script_name = 'treeshrink.py' 

def compose_build_distribution_name(build_type):
    return "treeshrink%s-v%s-%s" % (build_type, treeshrink.PROGRAM_VERSION, datetime.now().strftime("%Y%b%d"))

param = {
    'name': treeshrink.PROGRAM_NAME,
    'version': treeshrink.PROGRAM_VERSION,
    'description': treeshrink.PROGRAM_DESCRIPTION,
    'author': treeshrink.PROGRAM_AUTHOR,
    'author_email': ['umai@ucsd.edu'],
    'url': treeshrink.PROGRAM_WEBSITE,
    'license': treeshrink.PROGRAM_LICENSE,
    'packages': find_packages(),
    'package_dir': {'treeshrink': 'treeshrink'},
    'test_suite': "treeshrink.test",
    'include_package_data': True,
    'install_requires': ['dendropy==4.3.0'],
    'scripts' : [script_name],
    'zip_safe': True,
    'keywords': 'Phylogenetics Evolution Biology',
    'long_description': treeshrink.PROGRAM_DESCRIPTION,
    'classifiers': ["Environment :: Console",
                    "Intended Audience :: Developers",
                    "Intended Audience :: Science/Research",
                    "Natural Language :: English",
                   "Operating System :: OS Independent",                                                     "Programming Language :: Python",                                                         "Topic :: Scientific/Engineering :: Bio-Informatics",                                                                                                                                                               ],
    }

setup(**param)    
