import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from landaupy import langauss
from scipy.stats import median_abs_deviation, norm
from scipy.optimize import curve_fit

import ROOT as root

def binned_fit_langauss(samples, bins, min_x_val, max_x_val, nan='remove'):
  if nan == 'remove':
    samples = samples[~np.isnan(samples)]

  hist, bin_edges = np.histogram(samples, bins, range=(min_x_val,max_x_val), density=True)
  bin_centres = bin_edges[:-1] + np.diff(bin_edges) / 2
  hist = np.insert(hist, 0, sum(samples < bin_edges[0]))
  bin_centres = np.insert(bin_centres, 0, bin_centres[0] - np.diff(bin_edges)[0])
  hist = np.append(hist, sum(samples > bin_edges[-1]))
  bin_centres = np.append(bin_centres, bin_centres[-1] + np.diff(bin_edges)[0])

  hist = hist[1:-1]
  bin_centres = bin_centres[1:-1]
  landau_x_mpv_ansatz = bin_centres[np.argmax(hist)]
  landau_xi_ansatz = median_abs_deviation(samples) / 5
  gauss_sigma_ansatz = landau_xi_ansatz / 10

  popt, pcov = curve_fit(
    lambda x, mpv, xi, sigma: langauss.pdf(x, mpv, xi, sigma),
    xdata=bin_centres,
    ydata=hist,
    p0=[landau_x_mpv_ansatz, landau_xi_ansatz, gauss_sigma_ansatz],
  )
  return popt, pcov, hist, bin_centres

def gaussian(x, A, mu, sigma):
  return A * norm.pdf(x, mu, sigma)

def binned_fit_gaussian(samples, nBins, nan='remove'):
  if nan == 'remove':
    samples = samples[~np.isnan(samples)]
    samples = samples[~np.isinf(samples)]

  mu_ansatz = np.mean(samples)
  sigma_ansatz = np.std(samples)
  print(f"NUMBER OF SAMPLES: {len(samples)}")
  if mu_ansatz < 0:
    hist, bin_edges = np.histogram(samples, bins=nBins, range=(min(0, mu_ansatz - 2*sigma_ansatz), max(0, mu_ansatz + 2*sigma_ansatz)), density=True)
  else:
    print(min(2*mu_ansatz, mu_ansatz + 2*sigma_ansatz))
    hist, bin_edges = np.histogram(samples, bins=nBins, range=(max(0, mu_ansatz - 2*sigma_ansatz), min(2*mu_ansatz, mu_ansatz + 2*sigma_ansatz)), density=True)
  bin_centres = bin_edges[:-1] + np.diff(bin_edges) / 2

  A_ansatz = np.trapz(hist, bin_centres)
  print(mu_ansatz)
  print(sigma_ansatz)

  popt, pcov = curve_fit(gaussian, bin_centres, hist, p0=[A_ansatz, mu_ansatz, sigma_ansatz])
  A, mu, sigma = popt
  perr = np.sqrt(np.diag(pcov))

  y_fit = gaussian(bin_centres, A, mu, sigma)
  residuals = hist - y_fit

  sse = np.sum(residuals**2)
  sigma = np.sqrt(hist)
  sigma[sigma == 0] = 1
  chi2 = np.sum((residuals / sigma)**2)
  dof = len(hist) - len(popt)
  rchi2 = chi2 / dof

  popt_with_error = np.append(popt, perr[2])

  return popt_with_error, sse, rchi2

def get_fit_results_TR(df_of_results, mcp_specs):

  arr_of_sig_total_30 = [np.sqrt(sig_fit**2 + (1/600)) for sig_fit in df_of_results['Sigma']]
  arr_of_unc_total_30 = np.array(df_of_results['Uncertainty']) * np.array(df_of_results['Sigma']) / np.array(arr_of_sig_total_30)

  if (mcp_specs is not None) & (mcp_specs != (0, 0)):
    sig_dut_values_30 = []
    sig_dut_errors_30 = []
    mcp_tr = mcp_specs[0]/1000
    mcp_tr_err = mcp_specs[1]/1000
    print(f"[BETA ANALYSIS]: [TIME RESOLUTION] Calculating time resolution for DUT, assuming MCP time resolution {1000*mcp_tr} +/- {1000*mcp_tr_err} ps")
    for ch_ind, ch_val in enumerate(df_of_results['Sigma']):
      sig30 = np.sqrt(ch_val**2 - mcp_tr**2)
      sig30err = np.sqrt((ch_val*df_of_results['Uncertainty'][ch_ind])**2 + (mcp_tr*mcp_tr_err)**2)/sig30
      sig_dut_values_30.append((1000*sig30).round(1))
      sig_dut_errors_30.append((1000*sig30err).round(1))
    df_of_results['Resolution @ 30%'] = sig_dut_values_30
    df_of_results['Uncertainty @ 30%'] = sig_dut_errors_30

  else:
    print(f"[BETA ANALYSIS]: [TIME RESOLUTION] Calculating time resolution for a dual-plane setup, with unknown time resolutions")
    sig_30_p1 = np.sqrt(0.5*(df_of_results['Sigma'][0]**2 + df_of_results['Sigma'][2]**2 - df_of_results['Sigma'][1]**2))
    sig_30_p2 = np.sqrt(0.5*(df_of_results['Sigma'][1]**2 + df_of_results['Sigma'][2]**2 - df_of_results['Sigma'][0]**2))
    sig_30_p3 = np.sqrt(0.5*(df_of_results['Sigma'][0]**2 + df_of_results['Sigma'][1]**2 - df_of_results['Sigma'][2]**2))
    sig_30_err = np.sqrt(np.sum((np.array(df_of_results['Uncertainty'])*np.array(df_of_results['Sigma']))**2))
    sig_30_err1 = sig_30_err / sig_30_p1
    sig_30_err2 = sig_30_err / sig_30_p2
    sig_30_err3 = sig_30_err / sig_30_p3
    sig_dut_values_30 = [(1000*sig_30_p1).round(1),(1000*sig_30_p2).round(1),(1000*sig_30_p3).round(1)]
    sig_dut_errors_30 = [(1000*sig_30_err1).round(1),(1000*sig_30_err2).round(1),(1000*sig_30_err3).round(1)]
    df_of_results['Resolution @ 30%'] = sig_dut_values_30
    df_of_results['Uncertainty @ 30%'] = sig_dut_errors_30

  return df_of_results
