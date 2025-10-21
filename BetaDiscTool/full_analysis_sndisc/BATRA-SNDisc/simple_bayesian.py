import numpy as np
import sys

def gaussian(X, mu, sigma):
  return np.exp((-(X-mu)**2)/(2*sigma**2))/(np.sqrt(2*np.pi)*sigma)

def landau_approx(X, mu, sigma, kappa):
  return np.exp((-(X-mu))/(kappa) - np.exp((-(X-mu))/(kappa)))/(np.sqrt(2*np.pi)*sigma)

if len(sys.argv) != 3:
  print("Usage: python simple_bayesian.py <P_sig> <X_val>")
  sys.exit(1)

#P_sig = 0.2
#X_val = 10
P_sig = float(sys.argv[1])
X_val = float(sys.argv[2])
P_noise = 1 - P_sig

mu_s = 15
sigma_s = 1
kappa_s = 3

mu_n = 5
sigma_n = 5

P_X_given_s = landau_approx(X_val, mu_s, sigma_s, kappa_s)
P_X_given_n = gaussian(X_val, mu_n, sigma_n)
P_X = P_X_given_s * P_sig + P_X_given_n * P_noise

print(f"Amplitude {X_val} mV has event probability {round(P_X, 5)} given")
print(f"> signal: {round(P_X_given_s, 5)}")
print(f"> noise: {round(P_X_given_n, 5)}")
print(f"for 1 in {int(1/P_sig)} events being signal")
