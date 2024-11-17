#!/bin/bash

# This script is to be paired with the setup_custom_venv function found in notebook_utils.py

# Create a Python virtual environment
echo "Creating Python virtual environment..."
python3 -m venv custom_venv

# Activate the virtual environment
echo "Activating the virtual environment..."
source custom_venv/bin/activate

# # Install Tkinter for visualization
# echo "Installing Tkinter..."
# sudo apt-get install -y python3-tk

# Install required Python modules from the requirements file
echo "Installing required Python modules from requirements.txt..."
pip install -r requirements.txt

# Add the virtual environment as a kernel for Jupyter
echo "Adding virtual environment to Jupyter kernels..."
python -m ipykernel install --user --name=custom_venv --display-name "custom_venv"

# Provide user guidance for running the program
echo "Setup complete. To run the visualisation, use:"
echo "> source custom_venv/bin/activate"
echo "> python simple_NJ_tree.py"
