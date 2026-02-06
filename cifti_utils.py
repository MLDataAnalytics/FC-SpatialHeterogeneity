import os
import numpy as np
import nibabel as nib
from scipy import stats
from nilearn import image
import pandas as pd

class CiftiHandler:
    """Utility class for handling CIFTI files commonly used in the project."""
    
    def __init__(self, filepath=None):
        """Initialize CIFTI handler with optional file path."""
        self.filepath = filepath
        self.cifti_data = None
        if filepath:
            self.load(filepath)
    
    def load(self, filepath):
        """Load a CIFTI file."""
        try:
            self.cifti_data = nib.load(filepath)
            self.filepath = filepath
            return self.cifti_data
        except Exception as e:
            print(f"Error loading CIFTI file: {e}")
            return None
    
    def get_data(self):
        """Get the data array from CIFTI file."""
        if self.cifti_data is None:
            return None
        return self.cifti_data.get_fdata()
    
    def save(self, output_path):
        """Save CIFTI data to file."""
        if self.cifti_data is None:
            print("No data to save")
            return False
        try:
            nib.save(self.cifti_data, output_path)
            return True
        except Exception as e:
            print(f"Error saving CIFTI file: {e}")
            return False
    
    def get_atlas_labels(self, atlas_data):
        """Extract unique atlas labels from atlas data."""
        labels = np.unique(atlas_data[atlas_data > 0])
        return labels
    
    def extract_regional_values(self, data, atlas_data, labels):
        """Extract regional values from data using atlas parcellation."""
        regional_values = np.zeros(len(labels))
        for idx, label in enumerate(labels):
            mask = atlas_data == label
            if np.any(mask):
                regional_values[idx] = np.mean(data[mask])
        return regional_values


class StatisticalTests:
    """Utility class for statistical tests commonly used in the project."""
    
    @staticmethod
    def partial_correlation(x, y, covariates):
        """
        Calculate partial correlation between x and y, removing effect of covariates.
        
        Parameters:
        -----------
        x : array-like
            First variable
        y : array-like
            Second variable
        covariates : array-like
            Covariates to partial out (can be 2D for multiple covariates)
        
        Returns:
        --------
        float : partial correlation coefficient
        """
        x = np.asarray(x).flatten()
        y = np.asarray(y).flatten()
        covariates = np.atleast_2d(np.asarray(covariates))
        
        if covariates.ndim == 1:
            covariates = covariates.reshape(-1, 1)
        elif covariates.shape[0] == 1:
            covariates = covariates.T
        
        # Residualize x
        x_residuals = x - np.linalg.lstsq(covariates, x, rcond=None)[1]
        # Residualize y
        y_residuals = y - np.linalg.lstsq(covariates, y, rcond=None)[1]
        
        # Calculate correlation of residuals
        correlation = np.corrcoef(x_residuals, y_residuals)[0, 1]
        return correlation
    
    @staticmethod
    def ranksum_test(x, y, alternative='two-sided'):
        """
        Perform Mann-Whitney U test (Wilcoxon rank-sum test).
        
        Parameters:
        -----------
        x, y : array-like
            Input arrays
        alternative : str
            'two-sided', 'greater', or 'less'
        
        Returns:
        --------
        tuple : (statistic, p_value)
        """
        statistic, p_value = stats.mannwhitneyu(x, y, alternative=alternative)
        return statistic, p_value
    
    @staticmethod
    def benjamini_hochberg_fdr(p_values, fdr_threshold=0.05):
        """
        Apply Benjamini-Hochberg FDR correction.
        
        Parameters:
        -----------
        p_values : array-like
            P-values to correct
        fdr_threshold : float
            FDR threshold
        
        Returns:
        --------
        array : boolean array indicating significant values
        """
        p_values = np.asarray(p_values)
        sorted_idx = np.argsort(p_values)
        sorted_p = p_values[sorted_idx]
        
        # Calculate BH critical values
        n = len(p_values)
        bh_critical = (np.arange(1, n + 1) / n) * fdr_threshold
        
        # Find largest index where p <= critical value
        reject = np.zeros(n, dtype=bool)
        for i in range(n - 1, -1, -1):
            if sorted_p[i] <= bh_critical[i]:
                reject[sorted_idx[:i + 1]] = True
                break
        
        return reject
    
    @staticmethod
    def spin_test(observed_corr, permuted_corrs, tail='two-sided'):
        """
        Calculate p-value from spin test permutation.
        
        Parameters:
        -----------
        observed_corr : float
            Observed correlation
        permuted_corrs : array-like
            Correlations from permutations
        tail : str
            'two-sided', 'greater', or 'less'
        
        Returns:
        --------
        float : p-value
        """
        permuted_corrs = np.asarray(permuted_corrs).flatten()
        
        if tail == 'greater':
            p_value = np.sum(permuted_corrs >= observed_corr) / len(permuted_corrs)
        elif tail == 'less':
            p_value = np.sum(permuted_corrs <= observed_corr) / len(permuted_corrs)
        else:  # two-sided
            p_value = np.sum(np.abs(permuted_corrs) >= np.abs(observed_corr)) / len(permuted_corrs)
        
        return p_value


class DataIO:
    """Utility class for data input/output operations."""
    
    @staticmethod
    def load_mat_file(filepath, key=None):
        """Load MATLAB .mat file."""
        import scipy.io as sio
        try:
            mat_data = sio.loadmat(filepath)
            if key:
                return mat_data.get(key)
            return mat_data
        except Exception as e:
            print(f"Error loading MAT file: {e}")
            return None
    
    @staticmethod
    def load_csv(filepath, **kwargs):
        """Load CSV file using pandas."""
        try:
            return pd.read_csv(filepath, **kwargs)
        except Exception as e:
            print(f"Error loading CSV file: {e}")
            return None
    
    @staticmethod
    def save_csv(data, filepath, **kwargs):
        """Save data to CSV file."""
        try:
            if isinstance(data, pd.DataFrame):
                data.to_csv(filepath, **kwargs)
            else:
                np.savetxt(filepath, data, delimiter=',')
            return True
        except Exception as e:
            print(f"Error saving CSV file: {e}")
            return False


class VisualizationHelpers:
    """Helper functions for visualization."""
    
    @staticmethod
    def get_colormap(n_colors, cmap='tab10'):
        """Get colormap for visualization."""
        import matplotlib.cm as cm
        cmap_obj = cm.get_cmap(cmap)
        colors = [cmap_obj(i) for i in np.linspace(0, 1, n_colors)]
        return colors
    
    @staticmethod
    def create_figure_layout(n_plots, figsize=None):
        """Create figure layout for subplots."""
        import matplotlib.pyplot as plt
        n_cols = int(np.ceil(np.sqrt(n_plots)))
        n_rows = int(np.ceil(n_plots / n_cols))
        if figsize is None:
            figsize = (n_cols * 4, n_rows * 3)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        return fig, axes.flatten() if n_plots > 1 else np.array([axes])
