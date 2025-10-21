import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, poisson, median_abs_deviation
from scipy.optimize import curve_fit, minimize
from scipy.special import erf
import pandas as pd
from landaupy import langauss

import ROOT as root
from ROOT import TF1
import argparse
import glob
import re
import os
import csv
import sys

def binned_fit_langauss(samples, bins, min_x_val, max_x_val, channel, nan='remove'):
  if nan == 'remove':
    samples = samples[~np.isnan(samples)]

  hist, bin_edges = np.histogram(samples, bins, range=(min_x_val,max_x_val), density=True)
  bin_centres = bin_edges[:-1] + np.diff(bin_edges) / 2
  hist = np.insert(hist, 0, sum(samples < bin_edges[0]))
  bin_centres = np.insert(bin_centres, 0, bin_centres[0] - np.diff(bin_edges)[0])
  print(hist)
  hist = np.append(hist, sum(samples > bin_edges[-1]))
  bin_centres = np.append(bin_centres, bin_centres[-1] + np.diff(bin_edges)[0])

  hist = hist[1:-1]
  bin_centres = bin_centres[1:-1]
  landau_x_mpv_guess = bin_centres[np.argmax(hist)]
  landau_xi_guess = median_abs_deviation(samples) / 5
  gauss_sigma_guess = landau_xi_guess / 10

  popt, pcov = curve_fit(
    lambda x, mpv, xi, sigma: langauss.pdf(x, mpv, xi, sigma),
    xdata=bin_centres,
    ydata=hist,
    p0=[landau_x_mpv_guess, landau_xi_guess, gauss_sigma_guess],
  )
  return popt, pcov, hist, bin_centres

# Function to generate an approximate Langaus distribution
def langaus_sample(size, landau_loc=5, landau_scale=1.2, gauss_scale=0.5):
  landau_samples = np.random.standard_cauchy(size) * landau_scale + landau_loc  # Approximate Landau
  gaussian_samples = norm.rvs(scale=gauss_scale, size=size)  # Gaussian spread
  return landau_samples + gaussian_samples  # Convolution

def langaus(x, A, mu, sigma, k):
  landau_part = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-(x - mu) / k) * (1 + erf((x - mu) / (np.sqrt(2) * sigma)))
  return A * landau_part

def gaussian(x, A, mu, sigma):
  return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

'''
# Generate synthetic data
num_events = 10000
signal_amplitudes = langaus_sample(size=num_events // 2)  # Signal (Langaus)
noise_amplitudes = norm.rvs(size=num_events // 2, loc=2, scale=0.5)  # Noise (Gaussian)
max_amplitudes = np.concatenate([signal_amplitudes, noise_amplitudes])

signal_times = norm.rvs(size=num_events // 2, loc=50, scale=5)  # Gaussian for time (signal)
noise_times = np.random.uniform(0, 100, size=num_events // 2)  # Uniform noise times
max_times = np.concatenate([signal_times, noise_times])
'''

# Instead of using synthetic data, use the real data from ROOT files
# In this case, a good signal-noise-separated file from the DESY TB
# Run5_W16_W17_RT_Ex1/stats_Sr_Run5_Ch1-180V_Ch2-130V_Ch3-2750V_trig230V.root

config = [[3, 4.7], [1, 4.7], [1, 4.7], [2, 1]]
file_list = ['../../../../../DESY_TB_Dec_24/data_TB/Run4_W5_RT_Ex1/stats_Sr_Run4_Ch1-125V_Ch2-125V_Ch3-2750V_trig230V.root']
file_array = []
tree_array = []

for pattern in file_list:
  root_files = glob.glob(pattern)
  for root_file in root_files:
    try:
      theFile = root.TFile(root_file)
      file_array.append(theFile)
      tree_array.append(theFile.Get("Analysis"))
    except Exception as e:
      print(f"Error reading {root_file}: {e}")

for ch_ind, ch_val in enumerate(config):
  pmax_list = []
  width_list = []
  tmax_list = []
  if ch_val[0] == 1:
    for entry in tree_array[0]:
      pmax_sig = entry.pmax[ch_ind]
      negpmax_sig = entry.negpmax[ch_ind]
      tmax_sig = entry.tmax[ch_ind]
      width_sig = entry.width[ch_ind][2]
      A, B, C, D, E = 0, 1000, -100, -20, 20
      if (pmax_sig < A) or (pmax_sig > B) or (negpmax_sig < C) or (tmax_sig < D) or (tmax_sig > E):
        #print(f"BAD EVENTS {A} {B} {C} {D} {E}")
        continue
      else:
        pmax_list.append(pmax_sig)
        width_list.append(width_sig)
        tmax_list.append(tmax_sig)
    break
  else:
    continue

pmax = np.array(pmax_list)
tmax = np.array(tmax_list)
width = np.array(width_list)

ansatz_cut = 30
true_pmax_sig = pmax[pmax >= ansatz_cut]
true_pmax_noise = pmax[pmax < ansatz_cut]
true_width_sig = width[np.where(pmax >= ansatz_cut)[0]]
true_width_noise = width[np.where(pmax < ansatz_cut)[0]]
'''
plt.hist(true_pmax_noise, bins=200, range=(0,200), edgecolor='red')
plt.hist(true_pmax_sig, bins=200, range=(0,200), edgecolor='blue', alpha=0.6)
plt.xlabel("PMAX")
plt.ylabel("Frequency")
plt.yscale('log')
plt.show()

plt.hist(true_width_noise, bins=200, range=(0,10), edgecolor='red')
plt.hist(true_width_sig, bins=200, range=(0,10), edgecolor='blue', alpha=0.6)
plt.xlabel("Width")
plt.ylabel("Frequency")
plt.yscale('log')
plt.show()

plt.hist(true_pmax_noise/true_width_noise, bins=200, range=(0,200), edgecolor='red')
plt.hist(true_pmax_sig/true_width_sig, bins=200, range=(0,200), edgecolor='blue', alpha=0.6)
#plt.hist(pmax/width, bins=200, range=(0,200), edgecolor='black')
plt.xlabel("PMAX / Width")
plt.ylabel("Frequency")
plt.yscale('log')
plt.show()

sys.exit()
'''

num_events = len(pmax)
sig_events = np.count_nonzero(pmax >= ansatz_cut)
labels = np.concatenate([np.ones(sig_events), np.zeros(num_events - sig_events)])  # 1 = signal, 0 = noise

png_files = glob.glob("*.png")
for file in png_files:
  os.remove(file)

# Convert to PyTorch tensors
max_amplitudes = torch.tensor(pmax, dtype=torch.float32)
max_times = torch.tensor(tmax, dtype=torch.float32)
max_pow = torch.tensor(pmax / width, dtype=torch.float32)
labels = torch.tensor(labels, dtype=torch.float32)

# Define a differentiable signal probability model
class SignalProbabilityModel(nn.Module):
  def __init__(self):
    super().__init__()
    self.amp_mu = nn.Parameter(torch.tensor(2.0*ansatz_cut))   # Langaus center
    self.amp_sigma = nn.Parameter(torch.tensor(3.0)) # Langaus width
    self.time_mu = nn.Parameter(torch.tensor(-0.5))  # Gaussian center (signal-like time)
    self.time_sigma = nn.Parameter(torch.tensor(0.5)) # Gaussian width
    self.pow_mu = nn.Parameter(torch.tensor(2.0*ansatz_cut))  # Langaus center
    self.pow_sigma = nn.Parameter(torch.tensor(3.0)) # Langaus width

  def forward(self, amplitudes, times, pow):

    amp_prob = torch.sigmoid((amplitudes - self.amp_mu) / self.amp_sigma)
    time_prob = torch.exp(-0.5 * ((times - self.time_mu) / self.time_sigma) ** 2)
    pow_prob = torch.sigmoid((pow - self.pow_mu) / self.pow_sigma)
    return amp_prob * pow_prob  # Combined probability

#arr_num_epochs = [450,500,550,600,650,700,750,800]
#arr_prob_threshold = np.round(np.arange(0.35,0.81,0.01), 2) #[0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8]

precision = 0.005
plot_every = precision*1
dec_p = len(str(precision)) - 2
arr_num_epochs = [10000]
arr_prob_threshold = np.round(np.arange(0.04,0.08,precision), dec_p)

data_list = []
for num_epochs in arr_num_epochs:
  model = SignalProbabilityModel()
  optimizer = optim.Adam(model.parameters(), lr=0.01)

  for epoch in range(num_epochs):
    optimizer.zero_grad()
    scores = model(max_amplitudes, max_times, max_pow)
    loss = -torch.mean(labels * torch.log(scores + 1e-6) + (1 - labels) * torch.log(1 - scores + 1e-6))  
    loss.backward()
    optimizer.step()
    if epoch % 2000 == 0:
      print(f"Epoch {epoch}, Loss: {loss.item()}")
      #for name, param in model.named_parameters():
      #  print(name, param.grad)

  scores = model(max_amplitudes, max_times, max_pow).detach().numpy()
  df = pd.DataFrame(scores)
  df.to_csv("PMAX_TMAX_WIDTH30.csv", index=False, header=False)

  arr_of_MPV = []
  arr_of_width = []
  arr_of_sigma = []
  arr_mpv_frac = []
  arr_qmax_frac = []
  arr_of_sse = []
  arr_of_rchi2 = []
  plt.figure(figsize=(10, 6))

  for prob_threshold in arr_prob_threshold:
    selected_events = scores > prob_threshold
    filtered_amplitudes = max_amplitudes[selected_events].numpy()

    data_var = filtered_amplitudes
    xLower = 0
    xUpper = 150
    nBins = 150
    ch_ind = 1

    histo, bins, _ = plt.hist(data_var, bins=nBins, range=(xLower, xUpper), color='blue', edgecolor='black', alpha=0.6, density=True)
    print(len(data_var))

    bin_centres = bins[:-1] + np.diff(bins) / 2

    popt, pcov, fitted_hist, bin_centres = binned_fit_langauss(data_var, nBins, xLower, xUpper, ch_ind)
    arr_of_MPV.append(popt[0])
    arr_of_width.append(popt[1])
    arr_of_sigma.append(popt[2])

    count_1p0mpv = sum(1 for value in data_var if value > popt[0])
    count_1p5mpv = sum(1 for value in data_var if value > 1.5*popt[0])
    arr_mpv_frac.append(count_1p5mpv/count_1p0mpv)

    max_bin_index = np.argmax(histo)
    max_bin_centre = bin_centres[max_bin_index]
    count_max_bin = sum(1 for value in data_var if value > max_bin_centre)
    count_1p5max_bin = sum(1 for value in data_var if value > 1.5 * max_bin_centre)
    arr_qmax_frac.append(count_1p5max_bin / count_max_bin if count_max_bin > 0 else 0)

    mpv, xi, sigma = popt
    mu_lang, sigma_lang, k_lang = popt
    y_fit = langauss.pdf(bin_centres, mpv, xi, sigma)

    residuals = histo - y_fit

    sse = np.sum(residuals**2)
    norm_sse = sse / len(data_var)
    arr_of_sse.append(sse)
    sigma = np.sqrt(histo)
    sigma[sigma == 0] = 1
    chi2 = np.sum((residuals / sigma) ** 2)

    N = len(histo)
    p = len(popt)
    nu = N - p
    chi2_red = chi2 / nu
    print("SSE")
    print(sse)
    print("chi2")
    print(chi2)
    arr_of_rchi2.append(chi2_red)
    rchi2 = chi2_red
    print("red chi2")
    print(rchi2)

    bin_heights, bin_edges = np.histogram(max_amplitudes.numpy(), bins=150, range=(0, 150), density=True)
    bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
    popt_gaussian, _ = curve_fit(gaussian, bin_centres, bin_heights, p0=[1-(sig_events / num_events), 0.5*ansatz_cut, 3.0])
    A_gauss, mu_gauss, sigma_gauss = popt_gaussian
    filtered_heights, filtered_edges = np.histogram(filtered_amplitudes, bins=bin_edges, density=True)
    filtered_centres = (filtered_edges[:-1] + filtered_edges[1:]) / 2
    
    signal_event_count = len(filtered_amplitudes)
    total_event_count = len(max_amplitudes)
    scale_factor = signal_event_count / total_event_count
    langaus_scaled = langauss.pdf(bin_centres, *popt) * scale_factor

    '''
    popt_langaus, _ = curve_fit(langaus, bin_centres, filtered_heights, p0=[sig_events / num_events, 2.5*ansatz_cut, 3.0, 3.0])
    signal_event_count = len(filtered_amplitudes)
    total_event_count = len(max_amplitudes)
    scale_factor = signal_event_count / total_event_count
    langaus_scaled = langaus(bin_centres, *popt_langaus) * scale_factor
    A_lang, mu_lang, sigma_lang, k_lang = popt_langaus
    
    #residuals = filtered_heights - langaus(bin_centres, *popt_langaus)
    zero_for_no_events = filtered_heights
    zero_for_no_events[zero_for_no_events < 1/total_event_count] = 0
    residuals =  zero_for_no_events - langaus(bin_centres, *popt_langaus)
    sse = np.sum(residuals ** 2)
    norm_sse = sse / len(zero_for_no_events)
    sigma_for_chi2 = np.sqrt(zero_for_no_events)
    sigma_for_chi2[sigma_for_chi2 == 0] = 1
    chi2 = np.sum( (residuals / sigma_for_chi2) ** 2)
    rchi2 = chi2 / (len(bin_centres) - len(popt_langaus))
    print("NUMBER EVENTS")
    print(len(filtered_amplitudes))
    print("SSE")
    print(sse)
    print("Chi2")
    print(rchi2)
    '''

    if np.isclose(prob_threshold / plot_every, np.round(prob_threshold / plot_every)) == True:
      print(f"Number of epochs: {num_epochs} | Probability threshold: {prob_threshold}")
      fig, axes = plt.subplots(1, 2, figsize=(12, 5))
      axes[0].hist(max_amplitudes.numpy(), bins=150, range=(0, 150), alpha=0.7, label="Signal + Noise", color='lightgray', density=True)
      gaus_label = f"Gauss (Noise) Fit\nμ = {A_gauss:.2f}, σ = {mu_gauss:.2f}, k = {sigma_gauss:.2f}"
      axes[0].plot(bin_centres, gaussian(bin_centres, *popt_gaussian), color='orange', linestyle='-', label=gaus_label, linewidth=2)
      langaus_label1 = f"Rescaled Langaus (Signal) Fit\nμ = {mu_lang:.2f}, σ = {sigma_lang:.2f}, k = {k_lang:.2f}"
      axes[0].plot(bin_centres, langaus_scaled, 'g-', label=langaus_label1, linewidth=2)
      axes[0].set_xlabel("Max Amplitude")
      axes[0].set_ylabel("Normalised Counts")
      axes[0].set_title("Raw Signal + Noise Distribution")
      axes[0].set_yscale('log')
      axes[0].set_ylim(1/num_events, 1)
      axes[0].legend()
      
      axes[1].hist(filtered_amplitudes, bins=150, range=(0, 150), alpha=0.4, label="Filtered Signal", color='g', density=True)
      langaus_label = f"Langaus (Signal) Fit\nμ = {mu_lang:.2f}, σ = {sigma_lang:.2f}, k = {k_lang:.2f}\nRed. $\chi^{2}$ = {rchi2:.2e}"
      axes[1].plot(filtered_centres, langauss.pdf(filtered_centres, *popt), 'k--', label=langaus_label, linewidth=2)
      axes[1].set_xlabel("Max Amplitude")
      axes[1].set_ylabel("Normalised Counts")
      axes[1].set_title("Filtered Signal Distribution")
      axes[1].legend()
      axes[1].annotate(str(format(prob_threshold, "."+str(dec_p)+"f")), xy=(0.6,0.5), xycoords='axes fraction', fontsize=25, fontweight='bold')
      plt.tight_layout()
      plt.savefig("pmax_"+str(num_epochs)+"_"+str(format(prob_threshold, "."+str(dec_p)+"f"))[2:]+".png",facecolor='w')

    modchi2 = rchi2*pow(len(filtered_amplitudes),-1.4)
    data_list.append([num_epochs, prob_threshold, len(filtered_amplitudes), mu_lang.round(2), sse.round(6), norm_sse.round(10), chi2.round(2), rchi2.round(4), modchi2])

column_headings = ["Number of epochs","Probability threshold","Signal event count","MPV amplitude","SSE score","SSE / No. events","Chi2 value","Red. Chi2 value","Mod. Chi2 value"]
df = pd.DataFrame(data_list, columns=column_headings)
df.to_csv("disc_analysis.csv", index=False)

colours = ["r","orange","yellow","lime","green","blue","purple","magenta"]
markers = ["o","v","s","^","D","p","d","h"]

fig, axes = plt.subplots(1, 5, figsize=(25, 5))

for i, y_column in enumerate(['Signal event count', 'MPV amplitude', 'SSE score', 'Red. Chi2 value', 'Mod. Chi2 value']):
  ax = axes[i]
  for j, label in enumerate(df['Number of epochs'].unique()):
    subset = df[df['Number of epochs'] == label]
    ax.plot(subset['Probability threshold'], subset[y_column], label=f'{label} ({y_column})', color=colours[j], marker=markers[j], markersize=8, markeredgecolor='black', markeredgewidth=1)

  ax.set_xlabel('Probability threshold')
  ax.set_ylabel(y_column)
  if i >= 2:
    ax.set_yscale("log")
  ax.legend()

plt.tight_layout()
plt.savefig("pdperformanceplots.png",facecolor='w')

# Prefit selection
sse_thres = df['SSE score'].mean()
chi2_thres = 10E4
df = df.loc[df['SSE score'] <= sse_thres]
df = df.loc[df['Red. Chi2 value'] <= chi2_thres]

coeffs_quad = np.polyfit(df['Signal event count'], df['SSE score'], 2)
poly_eq = np.poly1d(coeffs_quad)
x_fit_quad = np.linspace(min(df['Signal event count']), max(df['Signal event count']), 100)
y_fit_quad = poly_eq(x_fit_quad)

def power_law(x, A):
  return x ** A

#A_opt, _ = curve_fit(power_law, df['Signal event count'], df['Red. Chi2 value'])
#x_fit_pow = np.linspace(min(df['Signal event count']), max(df['Signal event count']), 100)
#y_fit_pow = power_law(x_fit_pow, A_opt[0])

'''
# Post-fit selection
sse_thres = 1
chi2_thres = 10E4
df = df.loc[df['SSE score'] <= sse_thres]
df = df.loc[df['Red. Chi2 value'] <= chi2_thres]
'''

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for i, y_column in enumerate(['SSE score', 'Red. Chi2 value', 'Mod. Chi2 value']):
  ax = axes[i]
  ax.scatter(df['Signal event count'], df[y_column], c='g', s=1, marker='o')
  ax.set_xlabel('Signal event count')
  ax.set_ylabel(y_column)

  if i == 0:
    ax.plot(x_fit_quad, y_fit_quad, color='k', linestyle='--', label=f"Quadratic Fit: {coeffs_quad[0]:.2e}x² + {coeffs_quad[1]:.2e}x + {coeffs_quad[2]:.2e}")
    ax.legend()
  #if i == 1:
    #ax.plot(x_fit_pow, y_fit_pow, color='k', linestyle='--', label=f"Quadratic Fit: x^{A_opt[0]:.4f}")
    #ax.legend()

plt.tight_layout()
plt.savefig("pdscorevsize_linear.png",facecolor='w')

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for i, y_column in enumerate(['SSE score', 'Red. Chi2 value', 'Mod. Chi2 value']):
  ax = axes[i]
  ax.scatter(df['Signal event count'], df[y_column], c='g', s=1, marker='o')
  ax.set_xlabel('Signal event count')
  ax.set_ylabel(y_column)
  if i == 0: 
    ax.plot(x_fit_quad, y_fit_quad, color='k', linestyle='--', label=f"Quadratic Fit: {coeffs_quad[0]:.2e}x² + {coeffs_quad[1]:.2e}x + {coeffs_quad[2]:.2e}")
    ax.legend()
  #if i == 1: 
    #ax.plot(x_fit_pow, y_fit_pow, color='k', linestyle='--', label=f"Quadratic Fit: x^{A_opt[0]:.4f}")
    #ax.legend()
  ax.set_xscale("log")
  ax.set_yscale("log")

plt.tight_layout()
plt.savefig("pdscorevsize_log.png",facecolor='w')
