import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from landaupy import langauss

import ROOT as root
import argparse
import glob
import re
import os
import csv
import sys

from firstPass_ProcTools import getBias
from secondPass_ProcTools import binned_fit_langauss, gaussian
from diagnostics_ProcTools import perform_diagnostics
from class_SNDisc import SignalProbabilityModel, differential_programming_SigProbModel

def SNDisc_extract_signal(file, file_index, tree, channel_array, nBins, savename):

  png_files = glob.glob("NN_performance/*.png")
  for png_file in png_files:
    os.remove(png_file)
  max_pmax = nBins
  arr_signal_events = []
  amplitude_dfs = []

  for ch_ind, ch_val in enumerate(channel_array):
    pmax_list = []
    width_list = []
    tmax_list = []
    area_list = []
    risetime_list = []
    rms_list = []
    sensorType, AtQfactor, ansatz_pmax = ch_val
    if sensorType == 1:
      bias_of_channel = getBias(str(file), ch_ind)
      for entry in tree:
        pmax_sig = entry.pmax[ch_ind]
        width_sig = entry.width[ch_ind][4]
        #tmax_sig = entry.tmax[ch_ind]
        area_sig = entry.area_new[ch_ind]
        risetime_sig = entry.risetime[ch_ind]
        rms_sig = entry.rms[ch_ind]
        pmax_list.append(pmax_sig) # to fit Gaus noise and Langaus signal
        width_list.append(width_sig) # to fit Gaus noise and Langaus signal on pmax/width@30%
        #tmax_list.append(tmax_sig)
        area_list.append(area_sig)
        risetime_list.append(risetime_sig)
        rms_list.append(rms_sig)
    else:
      continue

    pmax = np.array(pmax_list)
    width = np.array(width_list)
    #tmax = np.array(tmax_list)
    area = np.array(area_list)
    risetime = np.array(risetime_list)
    rms = np.array(rms_list)

    num_epochs = 10000 # 10000
    precision = 0.0005 # 0.001
    learning_rate = 0.01
    plot_every = precision*1
    dec_p = len(str(precision)) - 2
    arr_prob_threshold = np.round(np.arange(0.01,0.12,precision), dec_p)

    scores, max_amplitudes = differential_programming_SigProbModel(num_epochs, pmax, width, area, risetime, rms, ansatz_pmax[file_index], learning_rate)
    df = pd.DataFrame(scores)
    df.to_csv(savename + "SNDisc_Ch"+str(ch_ind)+".csv", index=False, header=False)

    data_list = []
    plt.figure(figsize=(10, 6))

    scores[max_amplitudes.numpy() > max_pmax] = 0
    stopping_index_cond = len(arr_prob_threshold)
    print(f"Number of events in linear selection: {len(pmax[pmax > ansatz_pmax[file_index]])}")

    print(scores)

    for prob_threshold in arr_prob_threshold:
      selected_events = scores > prob_threshold
      filtered_amplitudes = max_amplitudes[selected_events].numpy()

      # Starting condition for probability threshold that doesn't cut away any signal events i.e. number of selected events is equal to the total
      if len(filtered_amplitudes) > 1.5*len(pmax[pmax > ansatz_pmax[file_index]]):
        continue

      # Stopping condition for probability threshold that there are fewer selected events by the NN than from the linear cut in PMAX
      if len(filtered_amplitudes) < 0.99*len(pmax[pmax > ansatz_pmax[file_index]]):
        stopping_index_cond = np.where(arr_prob_threshold == prob_threshold)[0]
        continue

      histo, bins, _ = plt.hist(filtered_amplitudes, bins=max_pmax, range=(0, max_pmax), color='blue', edgecolor='black', alpha=0.6, density=True)
      print(f"Number of filtered events: {len(filtered_amplitudes)}")

      bin_centres = bins[:-1] + np.diff(bins) / 2
      popt, pcov, fitted_hist, bin_centres = binned_fit_langauss(filtered_amplitudes, nBins, 0, max_pmax, ch_ind)

      mpv, xi, sigma = popt
      mu_lang, sigma_lang, k_lang = popt
      y_fit = langauss.pdf(bin_centres, mpv, xi, sigma)

      residuals = histo - y_fit
      sse = np.sum(residuals**2)
      norm_sse = sse / len(filtered_amplitudes)
      sigma = np.sqrt(histo)
      sigma[sigma == 0] = 1
      chi2 = np.sum((residuals / sigma) ** 2)

      nu = len(histo) - len(popt)
      chi2_red = chi2 / nu
      rchi2 = chi2_red

      bin_heights, bin_edges = np.histogram(max_amplitudes.numpy(), bins=max_pmax, range=(0, max_pmax), density=True)
      bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2
      popt_gaussian, _ = curve_fit(gaussian, bin_centres, bin_heights, p0=[1-(np.count_nonzero(pmax >= ansatz_pmax[file_index]) / len(pmax)), 0.5*ansatz_pmax[file_index], 0.5*ansatz_pmax[file_index]])
      A_gauss, mu_gauss, sigma_gauss = popt_gaussian
      filtered_heights, filtered_edges = np.histogram(filtered_amplitudes, bins=bin_edges, density=True)
      filtered_centres = (filtered_edges[:-1] + filtered_edges[1:]) / 2

      if np.isclose(prob_threshold / plot_every, np.round(prob_threshold / plot_every)) == True:
        signal_event_count = len(filtered_amplitudes)
        total_event_count = len(max_amplitudes)
        scale_factor = signal_event_count / total_event_count
        langaus_scaled = langauss.pdf(bin_centres, *popt) * scale_factor

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        axes[0].hist(max_amplitudes.numpy(), bins=max_pmax, range=(0, max_pmax), alpha=0.7, label="Signal + Noise", color='lightgray', density=True)
        gaus_label = f"Gauss (Noise) Fit\nμ = {A_gauss:.2f}, σ = {mu_gauss:.2f}, k = {sigma_gauss:.2f}"
        axes[0].plot(bin_centres, gaussian(bin_centres, *popt_gaussian), color='orange', linestyle='-', label=gaus_label, linewidth=2)
        langaus_label1 = f"Rescaled Langaus (Signal) Fit\nμ = {mu_lang:.2f}, k = {k_lang:.2f}, σ = {sigma_lang:.2f}"
        axes[0].plot(bin_centres, langaus_scaled, 'g-', label=langaus_label1, linewidth=2)
        axes[0].set_xlabel("Max Amplitude")
        axes[0].set_ylabel("Normalised Counts")
        axes[0].set_title("Raw Signal + Noise Distribution")
        axes[0].set_yscale('log')
        axes[0].set_ylim(1/len(pmax), 1)
        axes[0].legend()

        axes[1].hist(filtered_amplitudes, bins=max_pmax, range=(0, max_pmax), alpha=0.4, label="Filtered Signal", color='g', density=True)
        langaus_label = f"Langaus (Signal) Fit\nμ = {mu_lang:.2f}, σ = {sigma_lang:.2f}, k = {k_lang:.2f}\nRed. $\chi^{2}$ = {rchi2:.2e}"
        axes[1].plot(filtered_centres, langauss.pdf(filtered_centres, *popt), 'k--', label=langaus_label, linewidth=2)
        axes[1].set_xlabel("Max Amplitude")
        axes[1].set_ylabel("Normalised Counts")
        axes[1].set_title("Filtered Signal Distribution")
        axes[1].legend()
        axes[1].annotate(str(format(prob_threshold, "."+str(dec_p)+"f")), xy=(0.6,0.5), xycoords='axes fraction', fontsize=25, fontweight='bold')
        plt.tight_layout()
        plt.savefig("NN_performance/pmax_Ch"+str(ch_ind)+"_"+str(num_epochs)+"_"+str(format(prob_threshold, "."+str(dec_p)+"f"))[2:]+".png",facecolor='w')

      modchi2 = rchi2*pow(len(filtered_amplitudes),-1.4)
      num_ev_above_1p0 = len(filtered_amplitudes[filtered_amplitudes > mu_lang])
      num_ev_above_1p5 = len(filtered_amplitudes[filtered_amplitudes > 1.5*mu_lang])
      frac_ev_above_1p5 = num_ev_above_1p5 / num_ev_above_1p0
      data_list.append([ch_ind, bias_of_channel, num_epochs, prob_threshold, len(filtered_amplitudes), mu_lang.round(2), k_lang.round(3), sigma_lang.round(3), round(frac_ev_above_1p5,3), sse.round(6), norm_sse.round(10), chi2.round(4), rchi2.round(6), modchi2])

    column_headings = ["Channel","Bias","Number of epochs","Probability threshold","Signal event count","Amplitude MPV","Landau width","Gaussian sigma","Frac above 1p5 MPV","SSE score","SSE / No. events","Chi2","Red. Chi2","Mod. Chi2"]
    df = pd.DataFrame(data_list, columns=column_headings)
    change_in_number_events = df['Signal event count'].diff().dropna().to_numpy()
    #print("DIAGNOSTICS")
    #print(change_in_number_events)
    #print(arr_prob_threshold[1:stopping_index_cond+2])

    print("HERE IS SSE CALCS")
    print(df["SSE score"].dtype)
    print(df["SSE score"].unique())

    min_sse_prob = df.loc[df["SSE score"].idxmin(), "Probability threshold"]
    min_chi2_prob = df.loc[df["Red. Chi2"].idxmin(), "Probability threshold"]
    smallest_prob = min(min_sse_prob, min_chi2_prob)
    optimal_selected_events = scores > smallest_prob
    arr_signal_events.append(optimal_selected_events)

    np_ose = np.array(optimal_selected_events)
    print(f"Optimal number of filtered signal events: {len(np_ose[np_ose == True])}")
    print(f"Optimal cut at {smallest_prob}")

    amplitude_df_one_channel = df.loc[df["Probability threshold"] == smallest_prob]
    amplitude_df_one_channel = amplitude_df_one_channel[["Channel","Bias","Amplitude MPV","Landau width","Gaussian sigma","Frac above 1p5 MPV","SSE score","Red. Chi2"]]

    perform_diagnostics(df, ch_ind, savename)
    amplitude_dfs.append(amplitude_df_one_channel)

  amplitude_data = pd.concat(amplitude_dfs, ignore_index=True)
  return arr_signal_events, amplitude_data
