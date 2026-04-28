from setuptools import setup, find_packages
import treeshrink

param = {
    'name': treeshrink.PROGRAM_NAME,
    'version': treeshrink.PROGRAM_VERSION,
    'description': treeshrink.PROGRAM_DESCRIPTION,
    'author': treeshrink.PROGRAM_AUTHOR,
    'url': treeshrink.PROGRAM_WEBSITE,
    'license': treeshrink.PROGRAM_LICENSE,
    'packages': find_packages(),
    'include_package_data': True,
    'install_requires': [
        'treeswift',
        'numpy',
        'scipy',
    ],
    'python_requires': '>=3.8',
    'scripts' : ['run_treeshrink.py','decompose.py','make_gene_folder.py'],
    'keywords': 'Phylogenetics Evolution Biology',
    'long_description': """A Python implementation of the TreeShrink algorithm (Mai, Genome Biology, 2018)""",
    'long_description_content_type': 'text/plain',
    'classifiers': ["Environment :: Console",
                    "Intended Audience :: Developers",
                    "Intended Audience :: Science/Research",
                    "License :: OSI Approved :: GNU General Public License (GPL)",
                    "Natural Language :: English",
                    "Operating System :: OS Independent",
                    "Programming Language :: Python :: 3",
                    "Topic :: Scientific/Engineering :: Bio-Informatics",
                    ],
    }
    
setup(**param)
