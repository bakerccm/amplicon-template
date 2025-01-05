#!/bin/bash

# following instructions here:
# https://github.com/smirarab/sepp/tree/master/sepp-package

# starting from snakemake root
cd software

wget  "https://raw.github.com/smirarab/sepp-refs/master/gg/sepp-package.tar.bz"
tar xvfj sepp-package.tar.bz

cd sepp-package/sepp

# check python is available
# e.g. in base conda environment, this works fine
#    which python  # ~/miniconda3/bin/python
#    python --version  # Python 3.9.1

# setup script runs fast (~1 sec)
python setup.py config -c

echo "SEPP installation complete"
echo "You may now want to remove software/sepp-package.tar.bz"
