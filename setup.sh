#! /bin/bash

R CMD INSTALL dependencies/BMS_0.3.3.tar.gz

python setup.py install --user
