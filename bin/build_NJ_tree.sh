#!/bin/bash

# activate the virtual environment
source phylotree_venv/bin/activate

# install dependencies
pip install -r ../requirements/NJ_commandline.txt

# run the Python script and quickly visualise tree
python3 ../src/simple_NJ_tree.py
