from setuptools import setup, find_packages
import treeshrink
from os import walk, listdir
from os.path import join,normpath,isfile

def recursive_list_dir(path):
    listing=[]
    for x in walk(path):
        if isfile(x[0]):
            listing.append(x[0].split(path+'/')[1])
        for y in listdir(x[0]):
            z = normpath(join(x[0],y))
            if isfile(z):
                listing.append(z.split(path+'/')[1])
    return listing

param = {
    'name': treeshrink.PROGRAM_NAME,
    'version': treeshrink.PROGRAM_VERSION,
    'description': treeshrink.PROGRAM_DESCRIPTION,
    'author': treeshrink.PROGRAM_AUTHOR,
    'url': treeshrink.PROGRAM_WEBSITE,
    'license': treeshrink.PROGRAM_LICENSE,
    'packages': find_packages()+['Rlib','R_scripts'],
    'package_data':{'':recursive_list_dir('Rlib')+recursive_list_dir('R_scripts')},
    'include_package_data': True,
    'scripts' : ['run_treeshrink.py','decompose.py','make_gene_folder.py'],
    'keywords': 'Phylogenetics Evolution Biology',
    'long_description': """A Python implementation of the TreeShrink algorithm (Mai, Genome Biology, 2018)""",
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
