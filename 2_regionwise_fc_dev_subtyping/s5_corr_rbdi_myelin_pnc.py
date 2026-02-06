import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

class MyelinCorrelationAnalysis:
    def __init__(self, data):
        self.data = data

    def calculate_correlation(self, var1, var2):
        correlation, _ = pearsonr(self.data[var1], self.data[var2])
        return correlation

    def plot_correlation(self, var1, var2):
        plt.scatter(self.data[var1], self.data[var2])
        plt.xlabel(var1)
        plt.ylabel(var2)
        plt.title(f'Correlation between {var1} and {var2}')
        plt.show()

# Example usage:
# df = pd.read_csv('myelin_data.csv')
# myelin_analysis = MyelinCorrelationAnalysis(df)
# correlation = myelin_analysis.calculate_correlation('variable_1', 'variable_2')
# myelin_analysis.plot_correlation('variable_1', 'variable_2')
