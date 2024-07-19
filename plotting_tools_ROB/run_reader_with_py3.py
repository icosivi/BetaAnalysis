import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import poisson
import ROOT as root
from ROOT import TF1
from scipy.special import gammaln
import math
import pandas as pd
import argparse
import glob
import re
import os
import csv

#dir_name = "../TB_SPS_June_CMS_sensors/Run0_TB_SPS_LFoundry-K1/"

var_dict = {"tmax":"t_{max} / ns" , "pmax":"p_max / mV" , "negpmax":"-p_max / mV", "charge":"Q / fC", "area_new":"Area / pWb" , "rms":"RMS / mV"}

'''
time res v bias
collected charge v bias
rms v bias

t1-t2 Cividec
t3-t4 Minicircuit
'''

def get_fit_results(arr_of_results_to_fit,arr_of_biases,decomp_sigma=False,simplified=False):
  arr_of_mean = []
  arr_of_sigma = []
  arr_of_ampl = []
  arr_of_chi2 = []
  arr_of_ndf = []
  arr_of_red_chi2 = []
  #arr_of_prob = []

  for fit_func in arr_of_results_to_fit:
    mean = fit_func.GetParameter(1)  # Mean of the gauss distribution
    sigma = fit_func.GetParameter(2) # Sigma of the gauss distribution
    amplitude = fit_func.GetParameter(0)  # Amplitude of the gauss distribution
    chi2 = fit_func.GetChisquare()  # Chi-squared value of the fit
    ndf = fit_func.GetNDF()  # Number of degrees of freedom
    #prob = fit_func.GetProb()  # Probability of the fit result
    arr_of_mean.append(mean)
    arr_of_sigma.append(sigma)
    arr_of_ampl.append(amplitude)
    arr_of_chi2.append(chi2)
    arr_of_ndf.append(ndf)
    arr_of_red_chi2.append(chi2/ndf)
    #arr_of_prob.append(prob)

  df_of_results = pd.DataFrame({
    "Bias": arr_of_biases,
    "Mean": arr_of_mean,
    "Sigma": arr_of_sigma,
    "Amplitude": arr_of_ampl,
    "Chi2": arr_of_chi2,
    "NDF": arr_of_ndf,
    "RChi2": arr_of_red_chi2
    #"Prob": arr_of_prob
  })

  if decomp_sigma:
    sig1 = np.sqrt(0.5*(arr_of_sigma[0]**2 + arr_of_sigma[2]**2 - arr_of_sigma[1]**2))
    sig2 = np.sqrt(0.5*(arr_of_sigma[0]**2 + arr_of_sigma[1]**2 - arr_of_sigma[2]**2))
    sig3 = np.sqrt(0.5*(arr_of_sigma[1]**2 + arr_of_sigma[2]**2 - arr_of_sigma[0]**2))
    df_of_results['Sigma_cpt'] = ["sigma_1","sigma_2","sigma_3"]
    df_of_results['Sigma_value'] = [sig1,sig2,sig3]

  if simplified:
    if decomp_sigma:
      return df_of_results[["Bias","Sigma_cpt","Sigma_value"]]
    else:
      return df_of_results[["Bias","Mean","Sigma"]]

  else:
    return df_of_results

def getBias(filename):
  pattern = r"_(\d{2,3}V)."
  match = re.search(pattern, str(filename))
  if match:
    return str(match.group(1))
  else:
    return None

def hist_tree_file_basics(tree,file,var,index,nBins,xLower,xUpper,biasVal,cut_cond,ch,toaThreshold):
  if var == "cfd["+str(ch)+"]["+str(toaThreshold)+"]-cfd["+str(ch+1)+"]["+str(toaThreshold)+"]" or var == "cfd[0]["+str(toaThreshold)+"]-cfd["+str(ch)+"]["+str(toaThreshold)+"]":
    thisHist = root.TH1F("hist"+biasVal, var+";tn-tn+1 / ns ;Events", nBins, xLower, xUpper)
    tree.Draw(var+">>hist"+biasVal,cut_cond)
  elif var == "charge":
    thisHist = root.TH1F("hist"+biasVal, "area_new;"+var_dict[var]+";Events", nBins, xLower, xUpper)
    tree.Draw("area_new["+str(ch)+"]/4.7>>hist"+biasVal,cut_cond)
  else:
    thisHist = root.TH1F("hist"+biasVal, var+";"+var_dict[var]+";Events", nBins, xLower, xUpper)
    tree.Draw(var+"["+str(ch)+"]>>hist"+biasVal,cut_cond)
  thisHist.SetLineWidth(2)
  thisHist.SetLineColor(index+1)
  return thisHist

def plot_fit_curves(xLower,xUpper,fit_type,hist_to_fit,index,biasVal):
  thisFit = TF1(fit_type+"_hist"+biasVal, fit_type, xLower, xUpper)
  hist_to_fit.Fit(thisFit, "Q")
  thisFit.SetLineWidth(3)
  thisFit.SetLineColor(index+1)
  #thisFit.SetLineStyle(2)
  return thisFit

def tMaxTest(files,trees,ch,total_number_channels):
  print("tmax test run")
  var = "tmax"
  nBins = 1000
  cut_cond = "event>100"
  arr_of_hists = []
  arr_of_biases = []

  xLower = -2
  xUpper = 2

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      print(files[i])
      print(bias)
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch,0)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j,0)
        arr_of_hists.append(thisHist)
        arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of tmax distribution", 800, 600)
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(0, max_y)
  arr_of_hists[0].SetTitle("Comparison of tmax distribution")
  arr_of_hists[0].Draw()
  if len(arr_of_hists) > 1:
    for hist_to_draw in arr_of_hists[1:]:
      hist_to_draw.Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(i), "l")

  legend.Draw()

  c1.SaveAs("tmax_test.png")

def pMaxTest(files,trees,ch,total_number_channels):
  print("pmax test run")
  var = "pmax"
  nBins = 200
  arr_of_hists = []
  arr_of_biases = []
  cut_cond = "event>100"

  xLower = 0
  xUpper = 400

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch,0)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j,0)
        arr_of_hists.append(thisHist)
        arr_of_biases.append(bias)


  c1 = root.TCanvas("c1", "Comparison of pmax distribution", 800, 600)
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(0, max_y)
  arr_of_hists[0].SetTitle("Comparison of pmax distribution")
  arr_of_hists[0].Draw()
  for hist_to_draw in arr_of_hists:
    hist_to_draw.Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(i), "l")

  legend.Draw()

  c1.SaveAs("pmax_test.png")

def negPMax(files,trees,ch,total_number_channels):
  print("negpmax test run")
  var = "negpmax"
  nBins = 200
  arr_of_hists = []
  arr_of_biases = []
  cut_cond = "event>100"

  xLower = -40
  xUpper = 20

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch,0)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j,0)
        arr_of_hists.append(thisHist)
        arr_of_biases.append(bias)


  c1 = root.TCanvas("c1", "Comparison of negpmax distribution", 800, 600)
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(0, max_y)
  arr_of_hists[0].SetTitle("Comparison of negpmax distribution")
  arr_of_hists[0].Draw()
  for hist_to_draw in arr_of_hists:
    hist_to_draw.Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(i), "l")

  legend.Draw()

  c1.SaveAs("negpmax_test.png")

def landauToArea(files,trees,ch,total_number_channels):
  print("Fitting Landau curve to area distribution")
  var = "area_new"
  nBins = 200
  xLower = 0
  xUpper = 300

  arr_of_hists = []
  arr_of_biases = []

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      cut_cond = "event>100 && negpmax["+str(ch)+"] > -30 && negpmax["+str(ch)+"] < 5 && pmax["+str(ch)+"] > 50 && pmax["+str(ch)+"] < 350"
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch,0)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        cut_cond = "event>100 && negpmax["+str(j)+"] > -30 && pmax["+str(j)+"] > 50 && pmax["+str(j)+"] < 350"
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j,0)
        arr_of_hists.append(thisHist)
        arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of area_new distribution", 800, 600)
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(0, max_y)
  arr_of_hists[0].SetTitle("Comparison of area_new distribution")
  arr_of_hists[0].Draw()
  for hist_to_draw in arr_of_hists[1:]:
    hist_to_draw.Draw("SAME")

  arr_of_fits = []
  for i in range(len(arr_of_hists)):
    thisFit = plot_fit_curves(xLower,xUpper,"landau",arr_of_hists[i],i,arr_of_biases[i])
    arr_of_fits.append(thisFit)
    thisFit.Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(i), "l")
  
  legend.Draw()
  c1.SaveAs("Landau_analysis.png")

  landau_results = get_fit_results(arr_of_fits,arr_of_biases)
  print(landau_results)

def landauToCharge(files, trees, ch, total_number_channels):
  print("Fitting Landau curve to charge distributions")
  var = "charge"
  nBins = 200
  xLower = 0
  xUpper = 60

  arr_of_hists = []
  arr_of_biases = []

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      cut_cond = "event>100 && negpmax["+str(ch)+"] > -30 && negpmax["+str(ch)+"] < 5 && pmax[" + str(ch) + "] > 50 && pmax[" + str(ch) + "] < 350"
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch,0)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        cut_cond = "event>100 && negpmax["+str(j)+"] > -30 && pmax["+str(j)+"] > 50 && pmax["+str(j)+"] < 350"
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j,0)
        arr_of_hists.append(thisHist)
        arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of charge distributions", 800, 600)

  arr_of_fits = []
  arr_of_divided_hists = []
  for i, hist in enumerate(arr_of_hists):
    hist_divided = hist.Clone(f"hist_divided_{i}")
    for bin in range(1, hist_divided.GetNbinsX() + 1):
      bin_content = hist_divided.GetBinContent(bin)
      hist_divided.SetBinContent(bin, bin_content / 4.7)

    arr_of_divided_hists.append(hist_divided)
    thisFit = plot_fit_curves(xLower, xUpper, "landau", hist_divided, i, arr_of_biases[i])
    arr_of_fits.append(thisFit)
    thisFit.Draw("SAME")

  max_y = max(hist.GetMaximum() for hist in arr_of_divided_hists) * 1.05
  arr_of_divided_hists[0].GetYaxis().SetRangeUser(0, max_y)
  arr_of_divided_hists[0].SetTitle("Comparison of charge distribution")
  arr_of_divided_hists[0].Draw()
  for hist_to_draw in arr_of_divided_hists[1:]:
    hist_to_draw.Draw("SAME")

  for hist_to_draw in arr_of_fits:
    hist_to_draw.Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(i), "l")

  legend.Draw()
  c1.SaveAs("Landau_charge_analysis.png")

  landau_results = get_fit_results(arr_of_fits,arr_of_biases)
  print(landau_results)

def rmsNoise(files,trees,ch,total_number_channels):

  print("Fitting Gaussian curve to RMS distribution")
  var = "rms"
  nBins = 200
  xLower = 0
  xUpper = 5

  arr_of_hists = []
  arr_of_biases = []

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      cut_cond = "event>100 && negpmax["+str(ch)+"] > -30 && negpmax["+str(ch)+"] < 5 && pmax[" + str(ch) + "] > 50 && pmax[" + str(ch) + "] < 350"
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch,0)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        cut_cond = "event>100 && negpmax["+str(j)+"] > -30 && pmax[" + str(j) + "] > 50 && pmax[" + str(j) + "] < 350"
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j,0)
        arr_of_hists.append(thisHist)
        arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of RMS distribution", 800, 600)
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(0, max_y)
  arr_of_hists[0].SetTitle("Comparison of RMS distribution")
  arr_of_hists[0].Draw()
  
  for hist_to_draw in arr_of_hists[1:]:
    hist_to_draw.Draw("SAME")

  arr_of_fits = []
  for i in range(len(arr_of_hists)):
    thisFit = plot_fit_curves(xLower,xUpper,"gaus",arr_of_hists[i],i,arr_of_biases[i])
    arr_of_fits.append(thisFit)
    thisFit.Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(i), "l")
  
  legend.Draw()
  c1.SaveAs("rms_Gauss_analysis.png")

  gauss_results = get_fit_results(arr_of_fits,arr_of_biases)
  print(gauss_results)

def timeRes(files,trees,ch,total_number_channels):
  print("Time resolution analysis")
  vars = []
  cut_conds = []
  toaThreshold = 1
  if ch == -1:
    for ch_it in range(total_number_channels-1):
      var = "cfd["+str(ch_it)+"]["+str(toaThreshold)+"]-cfd["+str(ch_it+1)+"]["+str(toaThreshold)+"]"
      cut_cond = "event>100 && negpmax["+str(ch_it)+"] > -30 && negpmax["+str(ch_it)+"] < 5 && pmax["+str(ch_it)+"] > 50 && pmax["+str(ch_it)+"] < 350 && negpmax["+str(ch_it+1)+"] > -30 && negpmax["+str(ch_it+1)+"] < 5 && pmax["+str(ch_it+1)+"] > 50 && pmax["+str(ch_it+1)+"] < 350 && cfd[0]["+str(toaThreshold)+"]>-2.2 && cfd[1]["+str(toaThreshold)+"]>-2.5 && cfd[2]["+str(toaThreshold)+"]>-2.0"
      vars.append(var)
      cut_conds.append(cut_cond)
    vars.append("cfd[0]["+str(toaThreshold)+"]-cfd["+str(total_number_channels-1)+"]["+str(toaThreshold)+"]")
    cut_conds.append("event>100 && negpmax["+str(total_number_channels-1)+"] > -30 && negpmax["+str(total_number_channels-1)+"] < 5 && pmax["+str(total_number_channels-1)+"] > 50 && pmax["+str(total_number_channels-1)+"] < 350 && negpmax[0] > -30 && negpmax[0] < 5 && pmax[0] > 50 && pmax[0] < 350 && cfd[0]["+str(toaThreshold)+"]>-2.2 && cfd[1]["+str(toaThreshold)+"]>-2.5 && cfd[2]["+str(toaThreshold)+"]>-2.0")
  else:
    var = "cfd["+str(ch)+"]["+str(toaThreshold)+"]-cfd["+str(ch+1)+"]["+str(toaThreshold)+"]"
    cut_cond = "event>100 && negpmax["+str(ch)+"] > -30 && negpmax["+str(ch)+"] < 5 && pmax["+str(ch)+"] > 50 && pmax["+str(ch)+"] < 350 && negpmax["+str(ch+1)+"] > -30 && negpmax["+str(ch+1)+"] < 5 && pmax["+str(ch+1)+"] > 50 && pmax["+str(ch+1)+"] < 350 && cfd[0]["+str(toaThreshold)+"]>-2.2 && cfd[1]["+str(toaThreshold)+"]>-2.5 && cfd[2]["+str(toaThreshold)+"]>-2.0"
    vars.append(var)
    cut_conds.append(cut_cond)

  nBins = 200
  xLower = -2
  xUpper = 2

  arr_of_hists = []
  arr_of_biases = []

  for i in range(len(trees)):
    for j in range(len(vars)):
      bias = getBias(files[i])
      thisHist = hist_tree_file_basics(trees[i],files[i],vars[j],j,nBins,xLower,xUpper,bias,cut_conds[j],j,toaThreshold)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of time resolution (CFD@"+str(toaThreshold+1)+"0%) distribution", 800, 600)
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(0, max_y)
  arr_of_hists[0].SetTitle("Comparison of time resolution (CFD@"+str(toaThreshold+1)+"0%) distribution")
  arr_of_hists[0].Draw("")
  for hist_to_draw in arr_of_hists[1:]:
    hist_to_draw.Draw("SAME")

  arr_of_fits = []
  for i in range(len(arr_of_hists)):
    thisFit = plot_fit_curves(xLower,xUpper,"gaus",arr_of_hists[i],i,arr_of_biases[i])
    thisFit.Draw("SAME")
    arr_of_fits.append(thisFit)

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(i), "l")

  legend.Draw()
  c1.SaveAs("time_res_analysis.png")

  gauss_results = get_fit_results(arr_of_fits,arr_of_biases,decomp_sigma=True)
  print(gauss_results)

def waveform(files,trees,ch,total_number_channels):
  print("Waveform analysis")
  t_data = [[] for _ in range(total_number_channels)]
  w_data = [[] for _ in range(total_number_channels)]
  event_limit = 10
  for i in range(len(trees)):
    tree = trees[i]
    if ch != -1:
      event_cycle = 0
      for entry in tree:
        if (event_cycle > event_limit): continue
        event_cycle += 1
        t = entry.t
        w = entry.w
        pmax = entry.pmax
        if 50 < pmax[ch] < 350:
          t_data[0].extend(t[ch])
          w_data[0].extend(w[ch])
    else:
      for ch_it in range(total_number_channels):
        event_cycle = 0
        for entry in tree:
          if (event_cycle > event_limit): continue
          event_cycle += 1
          t = entry.t
          w = entry.w
          pmax = entry.pmax
          if 50 < pmax[ch_it] < 350:
            t_data[ch_it].extend(t[ch_it])
            w_data[ch_it].extend(w[ch_it])
    
  t_data = [np.array(channel) for channel in t_data]
  w_data = [np.array(channel) for channel in w_data]

  fig, ax = plt.subplots(figsize=(10, 6))

  colors = plt.cm.viridis(np.linspace(0, 1, total_number_channels))

  for i in range(3):
    ax.scatter(t_data[i], w_data[i], label=f'Channel {i}', color=colors[i])

  ax.set_xlim(-1e-8,1e-8)
  ax.set_ylim(min(min(inner_array) for inner_array in w_data),max(max(inner_array) for inner_array in w_data))
  ax.set_xlabel('t / ns', fontsize=14)
  ax.set_ylabel('Amplitude / mV', fontsize=14)
  ax.set_title(f'Waveform comparison, first {event_limit} events per channel', fontsize=16)
  ax.legend(loc = "upper right")
  ax.grid(True)
  plt.tight_layout()
  plt.savefig("waveform_comparison.png",facecolor='w')

def makeCSV(files,trees,total_number_channels):
  print("CSV of information produced for plotting charge, RMS noise, and time resolution against bias")
  print("Charge collection")
  nBins = 200
  xLower = 0
  xUpper = 60
  arr_of_hists = []
  arr_of_biases = []
  for i in range(len(trees)):
    for j in range(total_number_channels):
      bias = getBias(files[i])
      cut_cond = "event>100 && negpmax["+str(j)+"] > -30 && pmax["+str(j)+"] > 50 && pmax["+str(j)+"] < 350"
      thisHist = hist_tree_file_basics(trees[i],files[i],"charge",j,nBins,xLower,xUpper,bias,cut_cond,j,0)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  arr_of_fits = []
  arr_of_divided_hists = []
  for i, hist in enumerate(arr_of_hists):
    hist_divided = hist.Clone(f"hist_divided_{i}")
    for bin in range(1, hist_divided.GetNbinsX() + 1):
      bin_content = hist_divided.GetBinContent(bin)
      hist_divided.SetBinContent(bin, bin_content / 4.7)
    arr_of_divided_hists.append(hist_divided)
    thisFit = plot_fit_curves(xLower, xUpper, "landau", hist_divided, i, arr_of_biases[i])
    arr_of_fits.append(thisFit)
  charge_col = get_fit_results(arr_of_fits,arr_of_biases,simplified=True)

  print("RMS noise")
  nBins = 200
  xLower = 0
  xUpper = 5
  arr_of_hists = []
  arr_of_biases = []
  for i in range(len(trees)):
    for j in range(total_number_channels):
      bias = getBias(files[i])
      cut_cond = "event>100 && negpmax["+str(j)+"] > -30 && pmax[" + str(j) + "] > 50 && pmax[" + str(j) + "] < 350"
      thisHist = hist_tree_file_basics(trees[i],files[i],"rms",j,nBins,xLower,xUpper,bias,cut_cond,j,0)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  arr_of_fits = []
  for i in range(len(arr_of_hists)):
    thisFit = plot_fit_curves(xLower,xUpper,"gaus",arr_of_hists[i],i,arr_of_biases[i])
    arr_of_fits.append(thisFit)
  rms_noise = get_fit_results(arr_of_fits,arr_of_biases,simplified=True)

  print("Time resolution")
  vars = []
  cut_conds = []
  toaThreshold = 1
  for ch_it in range(total_number_channels-1):
    var = "cfd["+str(ch_it)+"]["+str(toaThreshold)+"]-cfd["+str(ch_it+1)+"]["+str(toaThreshold)+"]"
    cut_cond = "event>100 && negpmax["+str(ch_it)+"] > -30 && negpmax["+str(ch_it)+"] < 5 && pmax["+str(ch_it)+"] > 50 && pmax["+str(ch_it)+"] < 350 && negpmax["+str(ch_it+1)+"] > -30 && negpmax["+str(ch_it+1)+"] < 5 && pmax["+str(ch_it+1)+"] > 50 && pmax["+str(ch_it+1)+"] < 350 && cfd[0]["+str(toaThreshold)+"]>-2.2 && cfd[1]["+str(toaThreshold)+"]>-2.5 && cfd[2]["+str(toaThreshold)+"]>-2.0"
    vars.append(var)
    cut_conds.append(cut_cond)
  vars.append("cfd[0]["+str(toaThreshold)+"]-cfd["+str(total_number_channels-1)+"]["+str(toaThreshold)+"]")
  cut_conds.append("event>100 && negpmax["+str(total_number_channels-1)+"] > -30 && negpmax["+str(total_number_channels-1)+"] < 5 && pmax["+str(total_number_channels-1)+"] > 50 && pmax["+str(total_number_channels-1)+"] < 350 && negpmax[0] > -30 && negpmax[0] < 5 && pmax[0] > 50 && pmax[0] < 350 && cfd[0]["+str(toaThreshold)+"]>-2.2 && cfd[1]["+str(toaThreshold)+"]>-2.5 && cfd[2]["+str(toaThreshold)+"]>-2.0")
  nBins = 200
  xLower = -2
  xUpper = 2
  arr_of_hists = []
  arr_of_biases = []
  for i in range(len(trees)):
    for j in range(len(vars)):
      bias = getBias(files[i])
      thisHist = hist_tree_file_basics(trees[i],files[i],vars[j],j,nBins,xLower,xUpper,bias,cut_conds[j],j,toaThreshold)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  arr_of_fits = []
  for i in range(len(arr_of_hists)):
    thisFit = plot_fit_curves(xLower,xUpper,"gaus",arr_of_hists[i],i,arr_of_biases[i])
    arr_of_fits.append(thisFit)
  time_res = get_fit_results(arr_of_fits,arr_of_biases,decomp_sigma=True,simplified=True)

  merged_CC_RMS = pd.merge(charge_col.reset_index(), rms_noise.reset_index(), on=['index', 'Bias'], how='inner')
  merged_df = pd.merge(merged_CC_RMS, time_res.reset_index(), on=['index', 'Bias'], how='inner')
  merged_df.rename(columns={'Mean_x': 'Charge', 'Sigma_x': 'err_Q', 'Mean_y': 'RMS', 'Sigma_y': 'err_rms'}, inplace=True)
  outcsvfile = 'fit_bias_data.csv'
  if os.path.exists(outcsvfile):
    existing_df = pd.read_csv(outcsvfile)
    existing_df_reset = existing_df.reset_index(drop=True)
    new_df_reset = merged_df.reset_index(drop=True)
    common_rows = new_df_reset.merge(existing_df_reset, on=['Bias'], how='inner')
    if len(common_rows) != 0:
      print("This bias point already exists in file '{outcsvfile}', check that there are no duplicate data points in the output file.")
    combined_df = pd.concat([existing_df, merged_df], ignore_index=True)
    combined_df.to_csv(outcsvfile, index=False)
    print(f"File '{outcsvfile}' exists, read and saved new data to it for use in bias plotter.")
  else:
    merged_df.to_csv(outcsvfile, index=False)
    print(f"File '{outcsvfile}' created for use in bias plotter.")


def main():
  parser = argparse.ArgumentParser(description='Read .root files into an array.')
  parser.add_argument('files', metavar='F', type=str, nargs='+',
                      help='List of .root files or wildcard pattern (*.root)')
  parser.add_argument('--ch', type=int, nargs='?', default=0,
                      help="For analysis of a specific channel, choose an integer in the range 1 to the total number of channels")
  parser.add_argument('--doTMAXTest', action='store_true', help='Enable the doTMAXTest option')
  parser.add_argument('--doPMAXTest', action='store_true', help='Enable the doPMAXTest option')
  parser.add_argument('--doNegPMAX', action='store_true', help='Enable the doNegPMAX option')
  parser.add_argument('--doLandauToArea', action='store_true', help='Enable the doLandauToArea option')
  parser.add_argument('--doChargeDist', action='store_true', help='Enable the doChargeDist option')
  parser.add_argument('--doRMSNoise', action='store_true', help='Enable the doRMSNoise option')
  parser.add_argument('--doTimeRes', action='store_true', help='Enable the doTimeRes option')
  parser.add_argument('--doWaveform', action='store_true', help='Enable the doWaveform option')
  parser.add_argument('--csvOut', action='store_true', help='Output important data for plots against bias')
  args = parser.parse_args()

  file_array = []
  tree_array = []

  for pattern in args.files:
    root_files = glob.glob(pattern)
    for root_file in root_files:
      try:
        theFile = root.TFile(root_file)
        file_array.append(theFile)
        tree_array.append(theFile.Get("Analysis"))
        print(f"Successfully read {root_file}")
      except Exception as e:
        print(f"Error reading {root_file}: {e}")

  if len(file_array) == 1:
    print(f"Single input ROOT file read.")
  else:
    print(f"Total {len(file_array)} input ROOT files read.")

  if args.ch == 0:
    tree_with_channels = theFile.Get("Analysis")
    branch_with_number_channels = tree_with_channels.GetBranch("t")
    outer_vectors = root.std.vector('std::vector<double>')()
    tree_with_channels.SetBranchAddress("t", outer_vectors)
    outer_vector_count = 0
    for i in range(tree_with_channels.GetEntries()):
      tree_with_channels.GetEntry(i)
      outer_vector_count += len(outer_vectors)
    total_number_channels = int(outer_vector_count / branch_with_number_channels.GetEntries())
    print(f"Total {total_number_channels} channels in the file.")
  else:
    print(f"Analysing specifically CH {args.ch}.")
    total_number_channels = 1


  if args.doTMAXTest: tMaxTest(file_array,tree_array,args.ch-1,total_number_channels)
  if args.doPMAXTest: pMaxTest(file_array,tree_array,args.ch-1,total_number_channels)
  if args.doNegPMAX: negPMax(file_array,tree_array,args.ch-1,total_number_channels)
  if args.doLandauToArea: landauToArea(file_array,tree_array,args.ch-1,total_number_channels)
  if args.doChargeDist: landauToCharge(file_array,tree_array,args.ch-1,total_number_channels)
  if args.doRMSNoise: rmsNoise(file_array,tree_array,args.ch-1,total_number_channels)
  if args.doTimeRes: timeRes(file_array,tree_array,args.ch-1,total_number_channels)
  if args.doWaveform: waveform(file_array,tree_array,args.ch-1,total_number_channels)
  if args.csvOut: makeCSV(file_array,tree_array,total_number_channels)

if __name__ == "__main__":
    main()
