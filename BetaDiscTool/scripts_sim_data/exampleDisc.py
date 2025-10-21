import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.special import erf
import pandas as pd

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

# Generate synthetic data
num_events = 10000
signal_amplitudes = langaus_sample(size=num_events // 2)  # Signal (Langaus)
noise_amplitudes = norm.rvs(size=num_events // 2, loc=2, scale=0.5)  # Noise (Gaussian)
max_amplitudes = np.concatenate([signal_amplitudes, noise_amplitudes])

signal_times = norm.rvs(size=num_events // 2, loc=50, scale=5)  # Gaussian for time (signal)
noise_times = np.random.uniform(0, 100, size=num_events // 2)  # Uniform noise times
max_times = np.concatenate([signal_times, noise_times])

labels = np.concatenate([np.ones(num_events // 2), np.zeros(num_events // 2)])  # 1 = signal, 0 = noise

print(max_amplitudes.sum())

# Convert to PyTorch tensors
max_amplitudes = torch.tensor(max_amplitudes, dtype=torch.float32)
max_times = torch.tensor(max_times, dtype=torch.float32)
labels = torch.tensor(labels, dtype=torch.float32)

# Define a differentiable signal probability model
class SignalProbabilityModel(nn.Module):
  def __init__(self):
    super().__init__()
    self.amp_mu = nn.Parameter(torch.tensor(5.0))   # Langaus center
    self.amp_sigma = nn.Parameter(torch.tensor(1.5)) # Langaus width
    self.time_mu = nn.Parameter(torch.tensor(50.0))  # Gaussian center (signal-like time)
    self.time_sigma = nn.Parameter(torch.tensor(5.0)) # Gaussian width

  def forward(self, amplitudes, times):
    # Probability of being in Langaus region
    amp_prob = torch.sigmoid((amplitudes - self.amp_mu) / self.amp_sigma)
    # Probability of being in Gaussian time region
    time_prob = torch.exp(-0.5 * ((times - self.time_mu) / self.time_sigma) ** 2)
    return amp_prob * time_prob  # Combined probability

#arr_num_epochs = [450,500,550,600,650,700,750,800]
#arr_prob_threshold = np.round(np.arange(0.35,0.81,0.01), 2) #[0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8]

arr_num_epochs = [500]
arr_prob_threshold = [0.5]

data_list = []
for num_epochs in arr_num_epochs:
  model = SignalProbabilityModel()
  optimizer = optim.Adam(model.parameters(), lr=0.01)

  for epoch in range(num_epochs):
    optimizer.zero_grad()
    scores = model(max_amplitudes, max_times)
    loss = -torch.mean(labels * torch.log(scores + 1e-6) + (1 - labels) * torch.log(1 - scores + 1e-6))  
    loss.backward()
    optimizer.step()
    if epoch % 50 == 0:
      print(f"Epoch {epoch}, Loss: {loss.item()}")
      for name, param in model.named_parameters():
        print(name, param.grad)

  scores = model(max_amplitudes, max_times).detach().numpy()

  for prob_threshold in arr_prob_threshold:
    selected_events = scores > prob_threshold
    filtered_amplitudes = max_amplitudes[selected_events].numpy()

    bin_heights, bin_edges = np.histogram(max_amplitudes.numpy(), bins=200, range=(0, 20), density=True)
    bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
    popt_gaussian, _ = curve_fit(gaussian, bin_centres, bin_heights, p0=[0.3, 2, 0.5])
    A_gauss, mu_gauss, sigma_gauss = popt_gaussian
    filtered_heights, filtered_edges = np.histogram(filtered_amplitudes, bins=bin_edges, density=True)
    filtered_centres = (filtered_edges[:-1] + filtered_edges[1:]) / 2
    popt_langaus, _ = curve_fit(langaus, bin_centres, filtered_heights, p0=[0.5, 5, 1.0, 1.0])
    signal_event_count = len(filtered_amplitudes)
    total_event_count = len(max_amplitudes)
    scale_factor = signal_event_count / total_event_count
    langaus_scaled = langaus(bin_centres, *popt_langaus) * scale_factor
    A_lang, mu_lang, sigma_lang, k_lang = popt_langaus
    
    sse = np.sum((filtered_heights - langaus(bin_centres, *popt_langaus)) ** 2)
    chi2 = np.sum((filtered_heights - langaus(bin_centres, *popt_langaus)) ** 2 / (langaus(bin_centres, *popt_langaus) + 1e-10))
    rchi2 = chi2 / (len(bin_centres) - len(popt_langaus))

    if np.isclose(prob_threshold / 0.05, np.round(prob_threshold / 0.05)) == True:
      print(f"Number of epochs: {num_epochs} | Probability threshold: {prob_threshold}")
      fig, axes = plt.subplots(1, 2, figsize=(12, 5))
      axes[0].hist(max_amplitudes.numpy(), bins=200, range=(0, 20), alpha=0.7, label="Signal + Noise", color='lightgray', density=True)
      gaus_label = f"Gauss (Noise) Fit\nμ = {A_gauss:.2f}, σ = {mu_gauss:.2f}, k = {sigma_gauss:.2f}"
      axes[0].plot(bin_centres, gaussian(bin_centres, *popt_gaussian), color='orange', linestyle='-', label=gaus_label, linewidth=2)
      langaus_label1 = f"Rescaled Langaus (Signal) Fit\nμ = {mu_lang:.2f}, σ = {sigma_lang:.2f}, k = {k_lang:.2f}"
      axes[0].plot(bin_centres, langaus_scaled, 'g-', label=langaus_label1, linewidth=2)
      axes[0].set_xlabel("Max Amplitude")
      axes[0].set_ylabel("Normalised Counts")
      axes[0].set_title("Raw Signal + Noise Distribution")
      axes[0].legend()
      
      axes[1].hist(filtered_amplitudes, bins=200, range=(0, 20), alpha=0.4, label="Filtered Signal", color='g', density=True)
      langaus_label = f"Langaus (Signal) Fit\nμ = {mu_lang:.2f}, σ = {sigma_lang:.2f}, k = {k_lang:.2f}\nRed. $\chi^{2}$ = {rchi2:.2e}"
      axes[1].plot(filtered_centres, langaus(filtered_centres, *popt_langaus), 'k--', label=langaus_label, linewidth=2)
      axes[1].set_xlabel("Max Amplitude")
      axes[1].set_ylabel("Normalised Counts")
      axes[1].set_title("Filtered Signal Distribution")
      axes[1].legend()
      plt.tight_layout()
      plt.savefig("pmax_"+str(num_epochs)+"_"+str(prob_threshold)[2:]+".png",facecolor='w')

    modchi2 = rchi2*pow(len(filtered_amplitudes),-1.4)
    data_list.append([num_epochs, prob_threshold, len(filtered_amplitudes), mu_lang.round(2), sse.round(3), chi2.round(1), rchi2.round(3), modchi2])

column_headings = ["Number of epochs","Probability threshold","Signal event count","MPV amplitude","SSE score","Chi2 value","Red. Chi2 value","Mod. Chi2 value"]
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
sse_thres = 1
chi2_thres = 10E4
df = df.loc[df['SSE score'] <= sse_thres]
df = df.loc[df['Red. Chi2 value'] <= chi2_thres]

coeffs_quad = np.polyfit(df['Signal event count'], df['SSE score'], 2)
poly_eq = np.poly1d(coeffs_quad)
x_fit_quad = np.linspace(min(df['Signal event count']), max(df['Signal event count']), 100)
y_fit_quad = poly_eq(x_fit_quad)

def power_law(x, A):
  return x ** A

A_opt, _ = curve_fit(power_law, df['Signal event count'], df['Red. Chi2 value'])
x_fit_pow = np.linspace(min(df['Signal event count']), max(df['Signal event count']), 100)
y_fit_pow = power_law(x_fit_pow, A_opt[0])

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
  print(len(df[y_column]))

  if i == 0:
    ax.plot(x_fit_quad, y_fit_quad, color='k', linestyle='--', label=f"Quadratic Fit: {coeffs_quad[0]:.2e}x² + {coeffs_quad[1]:.2e}x + {coeffs_quad[2]:.2e}")
    ax.legend()
  if i == 1:
    ax.plot(x_fit_pow, y_fit_pow, color='k', linestyle='--', label=f"Quadratic Fit: x^{A_opt[0]:.4f}")
    ax.legend()

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
  if i == 1: 
    ax.plot(x_fit_pow, y_fit_pow, color='k', linestyle='--', label=f"Quadratic Fit: x^{A_opt[0]:.4f}")
    ax.legend()
  ax.set_xscale("log")
  ax.set_yscale("log")

plt.tight_layout()
plt.savefig("pdscorevsize_log.png",facecolor='w')
