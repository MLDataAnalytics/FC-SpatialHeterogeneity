# Correlation Analysis between RBD index and S-A Axis

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
# Assuming data is in CSV format
# data = pd.read_csv('path_to_your_data.csv')

# Preprocessing steps for the data
# Example: Handling missing values, normalization etc.
# data.dropna(inplace=True)

# Calculate Correlation
# assuming 'RBD_index' and 'S_A_axis' are the column names in your dataframe
correlation_matrix = data[['RBD_index', 'S_A_axis']].corr()

# Display correlation matrix
print("Correlation Matrix:")
print(correlation_matrix)

# Visualization
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm')
plt.title('Correlation Heatmap between RBD index and S-A Axis')
plt.show()