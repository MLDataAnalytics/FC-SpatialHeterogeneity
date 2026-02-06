import os
import numpy as np
import scipy.io
import nibabel as nib
from nilearn import image

# Define paths
mat_file_path = 'path_to_your_mat_file.mat'
cifti_output_path = 'output_file_name.dscalar.nii'

# Load MAT files
mat_data = scipy.io.loadmat(mat_file_path)

# Extract prediction results
predictions = mat_data['predictions']  # update key based on actual data structure

# Function to compute metrics
def compute_metrics(pred, true):
    correlation = np.corrcoef(pred, true)[0, 1]
    mae = np.mean(np.abs(pred - true))
    return correlation, mae

# Example: Assuming true_values are defined
# correlation, mae = compute_metrics(predictions, true_values)

# CIFTI file handling example
cifti_data = nib.load(cifti_output_path)

# Perform regional analysis
# Assuming brain parcellations are defined

# Save metrics to CIFTI format
nib.save(cifti_data, cifti_output_path)

# More detailed regional performance analysis logic here...