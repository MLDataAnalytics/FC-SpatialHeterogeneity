import numpy as np
import pandas as pd

class PrincipalGradientCorrelation:
    def __init__(self, data):
        self.data = data

    def compute_gradients(self):
        self.gradients = np.gradient(self.data)
        return self.gradients

    def calculate_correlation(self, other):
        if self.gradients is None:
            raise ValueError("Gradients not computed. Call compute_gradients first.")
        correlation = np.corrcoef(self.gradients, other.gradients)
        return correlation

# Example usage:
if __name__ == '__main__':
    data1 = np.random.rand(100)
    data2 = np.random.rand(100)
    pgc1 = PrincipalGradientCorrelation(data1)
    pgc2 = PrincipalGradientCorrelation(data2)
    pgc1.compute_gradients()
    pgc2.compute_gradients()
    correlation_matrix = pgc1.calculate_correlation(pgc2)
    print(correlation_matrix)