import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

dataset = pd.read_csv('uudet.csv')

dataset2 = pd.read_csv('vanhat.csv')

X2 = dataset2['FileSize'].dropna().values
X1 = dataset['FileSize'].dropna().values

plt.figure(figsize=(10, 6))
plt.hist(X1, bins=30, alpha=0.5, label='Uudet', edgecolor='black', color='blue')
plt.hist(X2, bins=30, alpha=0.5, label='Vanhat', edgecolor='black', color='red')

plt.title('Distribution of FileSize')
plt.xlabel('FileSize')
plt.ylabel('Frequency')
plt.legend()
plt.grid(True)
plt.show()