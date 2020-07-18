#! /bin/bash

v=$1

sed -i -e "s/version:.*$/version: $v/g" treeshrink/meta.yaml
conda-build treeshrink/
conda convert --platform all //anaconda3/conda-bld/osx-64/treeshrink-$v*.tar.bz2 -o //anaconda3/conda-bld/
anaconda upload -u smirarab --force //anaconda3/conda-bld/*/treeshrink-$v*.tar.bz2
