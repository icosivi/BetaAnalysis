import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

file1 = "PMAX.csv"
file2 = "PMAX_TMAX.csv"
file3 = "PMAX_WIDTH30.csv"

data1 = np.loadtxt(file1, delimiter=",")
data2 = np.loadtxt(file2, delimiter=",")
data3 = np.loadtxt(file3, delimiter=",")

all_data = np.concatenate([data1, data2, data3])
x_min, x_max = all_data.min(), all_data.max()

fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True)

bins = np.linspace(0, 1, 1000)
axes[0].hist(data1, bins=bins, color='blue')
axes[1].hist(data2, bins=bins, color='green')
axes[2].hist(data3, bins=bins, color='red')

axes[0].set_title("PMAX[LANGAUS⊗GAUS︀]-only model")
axes[1].set_title("PMAX[LANGAUS⊗GAUS︀] + \nTMAX[GAUS⊗UNIFORM] model")
axes[2].set_title("PMAX[LANGAUS⊗GAUS︀] + \nPMAX/WIDTH@30%[LANGAUS⊗GAUS︀] model")

for ax in axes:
    ax.set_xlabel("Probability of signal event")
    ax.set_yscale('log')
    ax.set_ylim(0.5,2E4)

axes[0].set_ylabel("Frequency")
axes[0].set_xlim(0,0.3)
#axes[0].set_ylim(None, max(max(np.histogram(data1, bins=bins)[0]),max(np.histogram(data2, bins=bins)[0]),max(np.histogram(data3, bins=bins)[0])))
#axes[0].set_yscale('log')

plt.tight_layout()
plt.show()
