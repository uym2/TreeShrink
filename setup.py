from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages
import treeshrink
param = {
    'name': treeshrink.PROGRAM_NAME,
    'version': treeshrink.PROGRAM_VERSION,
    'description': treeshrink.PROGRAM_DESCRIPTION,
    'author': treeshrink.PROGRAM_AUTHOR,
    'author_email': ['treeshrink-users@googlegroups.com'],
    'url': treeshrink.PROGRAM_WEBSITE,
    'license': treeshrink.PROGRAM_LICENSE,
    'packages': find_packages(),
    'package_dir': {'treeshrink': 'treeshrink'},
    'include_package_data': True,
    'install_requires': ['dendropy>=4.00'],
    'scripts' : ['treeshrink.py'],
    #'zip_safe': True,
    'keywords': 'Phylogenetics Evolution Biology',
    'long_description': """A Python implementation of the Practical Alignment using SATe and Transitivity. 
    The package requires configuration to refer to third-party tools such as ClustalW2, MAFFT, MUCLE, OPAL, Prank, and RAxML, HMMMER,
    and the code is heavily based on SATe""",
    'classifiers': ["Environment :: Console",
                    "Intended Audience :: Developers",
                    "Intended Audience :: Science/Research",
                    "License :: OSI Approved :: GNU General Public License (GPL)",
                    "Natural Language :: English",
                    "Operating System :: OS Independent",
                    "Programming Language :: Python",
                    "Topic :: Scientific/Engineering :: Bio-Informatics",
                    ],
    }
    
setup(**param)
