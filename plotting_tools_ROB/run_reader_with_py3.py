import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import poisson
import ROOT as root
from ROOT import TF1
from scipy.special import gammaln
import math
from math import exp, sqrt, pi
import pandas as pd
import argparse
import glob
import re
import os
import csv
import math

#dir_name = "../TB_SPS_June_CMS_sensors/Run0_TB_SPS_LFoundry-K1/"

# DESY TB cut conditions
# Run0 S1 cut_cond = "event>-1 && pmax[1] > 8 && pmax[1] < 120 && negpmax[1] > -15"
# Run0 S2 cut_cond = "event>-1 && pmax[2] > 30 && pmax[2] < 120 && negpmax[2] > -15"

# Run1 S1 cut_cond = "event>-1 && pmax[1] > 8 && pmax[1] < 120 && negpmax[1] > -18"
# Run1 S2 cut_cond = "event>-1 && pmax[2] > 34 && pmax[2] < 120 && negpmax[2] > -12"

# Run2/3/7 S1 cut_cond = "event>-1 && pmax[1] > 20 && pmax[1] < 160 && negpmax[1] > -25"
# Run2/3/7 S2 cut_cond = "event>-1 && pmax[2] > 13 && pmax[2] < 120 && negpmax[2] > -15"

# Run4 S1 cut_cond = "event>-1 && pmax[1] > 15 && pmax[1] < 140 && negpmax[1] > -16"
# Run4 S2 cut_cond = "event>-1 && pmax[2] > 15 && pmax[2] < 120 && negpmax[2] > -10"

# Run5 S1 cut_cond = "event>-1 && pmax[1] > 15 && pmax[1] < 140 && negpmax[1] > -15"
# Run5 S2 cut_cond = "event>-1 && pmax[2] > 15 && pmax[2] < 140 && negpmax[2] > -18"

# Run6 S1 cut_cond = "event>-1 && pmax[1] > 14 && pmax[1] < 110 && negpmax[1] > -10"
# Run6 S2 cut_cond = "event>-1 && pmax[2] > 25 && pmax[2] < 120 && negpmax[2] > -15"

# Run8 S1 cut_cond = "event>-1 && pmax[1] > 22 && pmax[1] < 280 && negpmax[1] > -22"
# Run8 S2 cut_cond = "event>-1 && pmax[2] > 40 && pmax[2] < 280 && negpmax[2] > -22"

# Run9/10 S1 cut_cond = "event>-1 && pmax[1] >  && pmax[1] <  && negpmax[1] > -"
# Run9/10 S2 cut_cond = "event>-1 && pmax[2] >  && pmax[2] <  && negpmax[2] > -"

# Run11 S1 cut_cond = "event>-1 && pmax[1] > 28 && pmax[1] < 120 && negpmax[1] > -25"

# Run12 S1 cut_cond = "event>-1 && pmax[1] > 14 && pmax[1] < 160 && negpmax[1] > -10"
# Run12 S2 cut_cond = "event>-1 && pmax[2] > 15 && pmax[2] < 200 && negpmax[2] > -30"

# Run13 S1 cut_cond = "event>-1 && pmax[1] > 16 && pmax[1] < 120 && negpmax[1] > -25"
# Run13 S2 cut_cond = "event>-1 && pmax[2] > 35 && pmax[2] < 200 && negpmax[2] > -30"

# Run14 S1 cut_cond = "event>-1 && pmax[1] > 20 && pmax[1] < 200 && negpmax[1] > -25"
# Run14 S2 cut_cond = "event>-1 && pmax[2] > 20 && pmax[2] < 200 && negpmax[2] > -25"

#myGlobalCut = "event>-1"
myGlobalCut = "event>-1 && pmax[0] > 10 && pmax[0] < 150 && tmax[0] > 11 && tmax[0] < 14"

var_dict = {"tmax":"t_{max} / 10 ns" , "pmax":"p_max / mV" , "negpmax":"-p_max / mV", "charge":"Q / fC", "area_new":"Area / pWb" , "rms":"RMS / mV"}

'''
time res v bias
collected charge v bias
rms v bias

t1-t2 Cividec
t3-t4 Minicircuit
'''

def round_to_sig_figs(x, sig):
  if x == 0:
    return 0
  else:
    return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)

def get_fit_results(arr_of_results_to_fit,arr_of_biases,decomp_sigma=False,simplified=False,add_threshold_frac=None):
  arr_of_mean = []
  arr_of_sigma = []
  arr_of_ampl = []
  arr_more_sig = []
  #arr_of_chi2 = []
  #arr_of_ndf = []
  arr_of_red_chi2 = []
  #arr_of_prob = []

  for fit_func in arr_of_results_to_fit:
    if fit_func.GetNpar() > 3:
      gaus_sig = fit_func.GetParameter(3)
    else:
      gaus_sig = 0
    mean = fit_func.GetParameter(1)  # Mean of the gauss distribution
    sigma = fit_func.GetParameter(2) # Sigma of the gauss distribution
    amplitude = fit_func.GetParameter(0)  # Amplitude of the gauss distribution
    chi2 = fit_func.GetChisquare()  # Chi-squared value of the fit
    ndf = fit_func.GetNDF()  # Number of degrees of freedom
    #prob = fit_func.GetProb()  # Probability of the fit result
    arr_of_mean.append(round_to_sig_figs(mean,4))
    arr_of_sigma.append(round_to_sig_figs(sigma,4))
    arr_of_ampl.append(round_to_sig_figs(amplitude,3))
    arr_more_sig.append(round_to_sig_figs(gaus_sig,4))
    #arr_of_chi2.append(chi2)
    #arr_of_ndf.append(ndf)
    if ndf == 0:
      ndf = 1
    arr_of_red_chi2.append(round_to_sig_figs((chi2/ndf),3))
    #arr_of_prob.append(prob)

  if arr_more_sig[0] == 0:
    if add_threshold_frac is not None:
      print(add_threshold_frac)
      df_of_results = pd.DataFrame({
        "Bias": arr_of_biases[0],
        "Mean": arr_of_mean,
        "Landau width": arr_of_sigma,
        "Gaussian sigma": arr_more_sig,
        "Frac above 1p5": add_threshold_frac,
        "Amplitude": arr_of_ampl,
        "RChi2": arr_of_red_chi2
      })
    else:
      df_of_results = pd.DataFrame({
        "Bias": arr_of_biases[0],
        "Mean": arr_of_mean,
        "Sigma": arr_of_sigma,
        "Amplitude": arr_of_ampl,
        #"Chi2": arr_of_chi2,
        #"NDF": arr_of_ndf,
        "RChi2": arr_of_red_chi2
        #"Prob": arr_of_prob
      })
  else:
    df_of_results = pd.DataFrame({
      "Bias": arr_of_biases[0],
      "Mean": arr_of_mean,
      "Landau width": arr_of_sigma,
      "Gaussian sigma": arr_more_sig,
      "Amplitude": arr_of_ampl,
      "RChi2": arr_of_red_chi2
    })

  if decomp_sigma:
    if len(arr_of_sigma) < 3:
      sig1 = np.sqrt(arr_of_sigma[0]**2 - 0.010**2)
      df_of_results['Time_res'] = [sig1]
    else:
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
  print(filename)
  match = re.search(pattern, str(filename))
  if match:
    return str(match.group(1))
  else:
    pattern_hyp = r"-(\d{2,3}V)."
    match = re.search(pattern_hyp, str(filename))
    if match:
      return str(match.group(1))
    else:
      print("[GetBias] : BIAS NOT FOUND")
      return None

def hist_tree_file_basics(tree,file,var,index,nBins,xLower,xUpper,biasVal,cut_cond,ch,toaThreshold):
  if var == "cfd["+str(ch)+"]["+str(toaThreshold)+"]-cfd["+str(ch+1)+"]["+str(toaThreshold)+"]" or var == "cfd[0]["+str(toaThreshold)+"]-cfd["+str(ch)+"]["+str(toaThreshold)+"]":
    thisHist = root.TH1F("hist"+biasVal, var+";tn-tn+1 / ns ;Events", nBins, xLower, xUpper)
    tree.Draw(var+">>hist"+biasVal,cut_cond)
  elif var == "cfd["+str(ch)+"][2]-cfd["+str(ch)+"][4]":
    thisHist = root.TH1F("hist"+biasVal, var+";t(20%)-t(50%) / ns ;Events", nBins, xLower, xUpper)
    tree.Draw(var+">>hist"+biasVal,cut_cond)
  elif var == "cfd["+str(ch)+"]["+str(toaThreshold)+"]-cfd[1]["+str(toaThreshold)+"]":
    thisHist = root.TH1F("hist"+biasVal, var+";tn-tn+1 / ns ;Events", nBins, xLower, xUpper)
    tree.Draw(var+">>hist"+biasVal,cut_cond)
  elif var == "charge":
    thisHist = root.TH1F("hist"+biasVal, var+";"+var_dict[var]+"/4.7;Events", nBins, xLower, xUpper)
    nEntries = tree.GetEntries()
    for i in range(nEntries):
      tree.GetEntry(i)
      python_cut_cond = cut_cond.replace("&&", "and").replace("||", "or")
        
      if eval(python_cut_cond, {}, {"event": tree.event, "pmax": tree.pmax, "negpmax": tree.negpmax}):
        scaled_value = tree.area_new[ch] / 4.7
        thisHist.Fill(scaled_value)

      #if tree.event > -1 and 8 < tree.pmax[ch] < 120 and tree.negpmax[ch] > -18:
      #  scaled_value = tree.area_new[ch] / 4.7
      #  thisHist.Fill(scaled_value)
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
  nBins = 100
  cut_cond = "event>-1" # 100 for SPS TB 24
  arr_of_hists = []
  arr_of_biases = []

  xLower = -20
  xUpper = 20

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
  c1.SetLogy()
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(1, max_y)
  arr_of_hists[0].SetTitle("Comparison of tmax distribution")
  arr_of_hists[0].Draw()
  if len(arr_of_hists) > 1:
    for hist_to_draw in arr_of_hists[1:]:
      hist_to_draw.Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    if ch == -1: ch_num = i+1
    if ch != -1: ch_num = ch+1
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(ch_num), "l")

  legend.Draw()

  c1.SaveAs("tmax_test.png")

def pMaxTest(files,trees,ch,total_number_channels):
  print("pmax test run")
  var = "pmax"
  nBins = 160
  arr_of_hists = []
  arr_of_biases = []
  #cut_cond = "event>-1 && negpmax[0] > -15 && negpmax[0] < 5 && pmax[1] > 0 && pmax[1] < 120 && negpmax[1] > -15" # beta setup cut conditions
  cut_cond = "event>-1"

  xLower = 0
  xUpper = 160

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
  c1.SetLogy()
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(1, max_y)
  arr_of_hists[0].SetTitle("Comparison of pmax distribution")
  arr_of_hists[0].Draw()
  for hist_to_draw in arr_of_hists:
    hist_to_draw.Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    if ch == -1: ch_num = i+1
    if ch != -1: ch_num = ch+1
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(ch_num), "l")

  legend.Draw()

  #arr_of_fits = []
  #for i in range(len(arr_of_hists)):
  #  thisFit = plot_fit_curves(xLower,xUpper,"landau",arr_of_hists[i],i,arr_of_biases[i])
  #  arr_of_fits.append(thisFit)
  #  thisFit.Draw("SAME")

  c1.SaveAs("pmax_test.png")

  #pmax_landau = get_fit_results(arr_of_fits,arr_of_biases)
  #print(pmax_landau)

def negPMax(files,trees,ch,total_number_channels):
  print("negpmax test run")
  var = "negpmax"
  nBins = 200
  arr_of_hists = []
  arr_of_biases = []
  cut_cond = "event>-1"

  xLower = -35
  xUpper = 15

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
  c1.SetLogy()
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(1, max_y)
  arr_of_hists[0].SetTitle("Comparison of negpmax distribution")
  arr_of_hists[0].Draw()
  for hist_to_draw in arr_of_hists:
    hist_to_draw.Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    if ch == -1: ch_num = i+1
    if ch != -1: ch_num = ch+1
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(ch_num), "l")

  legend.Draw()

  c1.SaveAs("negpmax_test.png")

def landauToArea(files,trees,ch,total_number_channels):
  print("Fitting Landau curve to area distribution")
  var = "area_new"
  nBins = 80
  xLower = 0
  xUpper = 80

  arr_of_hists = []
  arr_of_biases = []

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      #cut_cond = "event>-1 && negpmax["+str(ch)+"] > -15 && negpmax["+str(ch)+"] < 5 && pmax["+str(ch)+"] > 6.0 && pmax["+str(ch)+"] < 120"
      cut_cond = myGlobalCut
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch,0)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        #cut_cond = "event>-1 && negpmax["+str(j)+"] > -15 && pmax["+str(j)+"] > 6.0 && pmax["+str(j)+"] < 120"
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j,0)
        arr_of_hists.append(thisHist)
        arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of area distribution", 800, 600)
  #c1.SetLogy()
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(1, max_y)
  arr_of_hists[0].SetTitle("Comparison of area distribution")
  arr_of_hists[0].Draw()

  for hist_to_draw in arr_of_hists[1:]:
    hist_to_draw.Draw("SAME")

  arr_of_fits = []
  for i in range(len(arr_of_hists)):
    print(arr_of_hists[i].GetMaximum())
    thisFit = plot_fit_curves(xLower, xUpper, "landau", arr_of_hists[i], i, arr_of_biases[i])
    arr_of_fits.append(thisFit)
    MPV = thisFit.GetParameter(1)
    thisFit.Draw("SAME")
    break

  total_events_above_threshold1p5 = 0
  total_events_above_threshold1p0 = 0
  for bin in range(1, arr_of_hists[0].GetNbinsX() + 1):
    bin_centre = arr_of_hists[0].GetBinCenter(bin)
    bin_content = arr_of_hists[0].GetBinContent(bin)
    if bin_centre > MPV:
      total_events_above_threshold1p0 += bin_content
    if bin_centre > 1.5*MPV:
      total_events_above_threshold1p5 += bin_content

  print(thisFit.GetMaximum(xLower,xUpper))
  #print(arr_of_hists[0].GetBinWidth(1))
  #print(thisFit.GetParameter(0) / arr_of_hists[0].GetBinWidth(1))

  arr_of_fits[0].Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    if ch == -1: ch_num = i+1
    if ch != -1: ch_num = ch+1
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(ch_num), "l")
  
  legend.Draw()
  c1.SaveAs("Landau_analysis.png")

  frac_thres = total_events_above_threshold1p5 / total_events_above_threshold1p0
  print(frac_thres)
  landau_results = get_fit_results(arr_of_fits,arr_of_biases,add_threshold_frac=frac_thres)
  print(landau_results)

def landauToCharge(files, trees, ch, total_number_channels):
  print("Fitting Landau curve to charge distributions")
  var = "charge"
  nBins = 60
  xLower = 0
  xUpper = 60

  arr_of_hists = []
  arr_of_biases = []

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      #cut_cond = "event>-1 && negpmax["+str(ch)+"] > -15 && negpmax["+str(ch)+"] < 5 && pmax[" + str(ch) + "] > 6.0 && pmax[" + str(ch) + "] < 120" # beta setup cut conds
      cut_cond = myGlobalCut
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch,0)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        #cut_cond = "event>-1 && negpmax["+str(j)+"] > -15 && pmax["+str(j)+"] > 6.0 && pmax["+str(j)+"] < 200 && pmax[1] > 0 && pmax[1] < 120 && negpmax[1] > -18" # beta setup cut conds
        cut_cond = "event>-1"
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j,0)
        arr_of_hists.append(thisHist)
        arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of charge distributions", 800, 600)
  c1.SetLogy()

  arr_of_fits = []
  for i, hist in enumerate(arr_of_hists):
    thisFit = plot_fit_curves(xLower, xUpper, "landau", hist, i, arr_of_biases[i])
    arr_of_fits.append(thisFit)
    MPV = thisFit.GetParameter(1)

  total_events_above_threshold1p5 = 0
  total_events_above_threshold1p0 = 0
  
  for bin in range(1, arr_of_hists[0].GetNbinsX() + 1):
    bin_centre = arr_of_hists[0].GetBinCenter(bin)
    bin_content = arr_of_hists[0].GetBinContent(bin)
    if bin_centre > MPV:
      total_events_above_threshold1p0 += bin_content
    if bin_centre > 1.5*MPV:
      total_events_above_threshold1p5 += bin_content
  
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(1, max_y)
  arr_of_hists[0].SetTitle("Comparison of charge distribution")
  arr_of_hists[0].Draw()
  for hist_to_draw in arr_of_hists[1:]:
    hist_to_draw.Draw("SAME")

  for thisFit in arr_of_fits:
    thisFit.Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    if ch == -1: ch_num = i+1
    if ch != -1: ch_num = ch+1
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(ch_num), "l")

  legend.Draw()
  c1.SaveAs("charge_analysis.png")

  frac_thres = total_events_above_threshold1p5 / total_events_above_threshold1p0
  landau_results = get_fit_results(arr_of_fits,arr_of_biases,add_threshold_frac=frac_thres)
  print(landau_results)

def rmsNoise(files,trees,ch,total_number_channels):

  print("Fitting Gaussian curve to RMS distribution")
  var = "rms"
  nBins = 100
  xLower = 0
  xUpper = 5

  arr_of_hists = []
  arr_of_biases = []

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      #cut_cond = "event>-1 && negpmax["+str(ch)+"] > -15 && negpmax["+str(ch)+"] < 5 && pmax[" + str(ch) + "] > 6.0 && pmax[" + str(ch) + "] < 120 && pmax[1] > 0 && pmax[1] < 120" # beta setup cut cond
      cut_cond = myGlobalCut
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch,0)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        #cut_cond = "event>-1 && negpmax["+str(j)+"] > -15 && pmax[" + str(j) + "] > 6.0 && pmax[" + str(j) + "] < 120 && pmax[1] > 0 && pmax[1] < 120 && negpmax[1] > -15" # beta setup cut cond
        cut_cond = "event>-1"
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j,0)
        arr_of_hists.append(thisHist)
        arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of RMS distribution", 800, 600)
  c1.SetLogy()
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(1, max_y)
  arr_of_hists[0].SetTitle("Comparison of RMS distribution")
  arr_of_hists[0].Draw()
  
  for hist_to_draw in arr_of_hists[1:]:
    hist_to_draw.Draw("SAME")

  arr_of_fits = []
  for i in range(len(arr_of_hists)):
    thisFit = plot_fit_curves(xLower,xUpper,"gaus",arr_of_hists[i],i,arr_of_biases[i])
    arr_of_fits.append(thisFit)
    thisFit.Draw("SAME")
    #break

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    if ch == -1: ch_num = i+1
    if ch != -1: ch_num = ch+1
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(ch_num), "l")
  
  legend.Draw()
  c1.SaveAs("rms_Gauss_analysis.png")

  gauss_results = get_fit_results(arr_of_fits,arr_of_biases)
  print(gauss_results)

def timeRes(files,trees,ch,total_number_channels):
  print("Time resolution analysis")
  vars = []
  cut_conds = []
  toaThreshold = 2
  # 1 = 20%, 4 = 50%
  # "wt" Draw option to ignore tail
  if total_number_channels == 2:
    var = "cfd[0]["+str(toaThreshold)+"]-cfd[2]["+str(toaThreshold)+"]"
    cut_cond = "negpmax[0] > -15 && negpmax[0] < 5 && pmax[0] > 6.0 && pmax[0] < 120 && pmax[1] > 20 && pmax[1] < 120 && negpmax[1] > -15 && cfd[1]["+str(toaThreshold)+"]<2 && cfd[1]["+str(toaThreshold)+"]>-2 && cfd[0]["+str(toaThreshold)+"]>10"
    vars.append(var)
    cut_conds.append(cut_cond)
  elif total_number_channels == 1:
    var = "cfd["+str(ch)+"]["+str(toaThreshold)+"]-cfd[1]["+str(toaThreshold)+"]"
    cut_cond = "negpmax[1] > -40 && pmax[1] < 400 && " + myGlobalCut
    vars.append(var)
    cut_conds.append(cut_cond)
  else:
    if ch == -1:
      #for ch_it in range(total_number_channels-1):
      for ch_it in range(3): # use this to force to look at first 3 channels, relevant for the 4-channel TRG-DUT-DUT-MCP test beam setup
        var = "cfd["+str(ch_it)+"]["+str(toaThreshold)+"]-cfd["+str(ch_it+1)+"]["+str(toaThreshold)+"]"
        #cut_cond = "event>-1 && negpmax["+str(ch_it)+"] > -15 && negpmax["+str(ch_it)+"] < 5 && pmax["+str(ch_it)+"] > 6.0 && pmax["+str(ch_it)+"] < 120 && negpmax["+str(ch_it+1)+"] > -15 && negpmax["+str(ch_it+1)+"] < 5 && pmax["+str(ch_it+1)+"] > 6.0 && pmax["+str(ch_it+1)+"] < 120 && cfd[0]["+str(toaThreshold)+"]>-2.2 && cfd[1]["+str(toaThreshold)+"]>-2.5 && cfd[2]["+str(toaThreshold)+"]>-2.0" # beta setup cut conds
        cut_cond = "event>-1"
        vars.append(var)
        cut_conds.append(cut_cond)
      vars.append("cfd[0]["+str(toaThreshold)+"]-cfd["+str(total_number_channels-1)+"]["+str(toaThreshold)+"]")
      cut_conds.append("event>-1 && negpmax["+str(total_number_channels-1)+"] > -15 && negpmax["+str(total_number_channels-1)+"] < 5 && pmax["+str(total_number_channels-1)+"] > 6.0 && pmax["+str(total_number_channels-1)+"] < 120 && negpmax[0] > -15 && negpmax[0] < 5 && pmax[0] > 6.0 && pmax[0] < 120 && cfd[0]["+str(toaThreshold)+"]>-2.2 && cfd[1]["+str(toaThreshold)+"]>-2.5 && cfd[2]["+str(toaThreshold)+"]>-2.0")
    else:
      var = "cfd["+str(ch)+"]["+str(toaThreshold)+"]-cfd["+str(ch+1)+"]["+str(toaThreshold)+"]"
      cut_cond = "event>-1 && negpmax["+str(ch)+"] > -15 && negpmax["+str(ch)+"] < 5 && pmax["+str(ch)+"] > 6.0 && pmax["+str(ch)+"] < 120 && negpmax["+str(ch+1)+"] > -15 && negpmax["+str(ch+1)+"] < 5 && pmax["+str(ch+1)+"] > 6.0 && pmax["+str(ch+1)+"] < 120 && cfd[0]["+str(toaThreshold)+"]>-2.2 && cfd[1]["+str(toaThreshold)+"]>-2.5 && cfd[2]["+str(toaThreshold)+"]>-2.0"
      vars.append(var)
      cut_conds.append(cut_cond)

  nBins = 80
  xLower = 8
  xUpper = 12

  arr_of_hists = []
  arr_of_biases = []

  print(vars)
  print(cut_conds)

  if total_number_channels >= 3: range_of_vars = 3
  else: range_of_vars = len(vars)

  for i in range(len(trees)):
    for j in range(range_of_vars):
      bias = getBias(files[i])
      if range_of_vars == 1:
        thisHist = hist_tree_file_basics(trees[i],files[i],vars[j],j,nBins,xLower,xUpper,bias,cut_conds[j],ch,toaThreshold)
      else:
        thisHist = hist_tree_file_basics(trees[i],files[i],vars[j],j,nBins,xLower,xUpper,bias,cut_conds[j],j,toaThreshold)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)

  print(len(arr_of_hists))
  print(len(arr_of_biases))

  c1 = root.TCanvas("c1", "Comparison of time resolution (CFD@"+str(toaThreshold+1)+"0%) distribution", 800, 600)
  c1.SetLogy()
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(1, max_y)
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
    if ch == -1: ch_num = i+1
    if ch != -1: ch_num = ch+1
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(ch_num), "l")

  legend.Draw()
  c1.SaveAs("time_res_analysis.png")

  gauss_results = get_fit_results(arr_of_fits,arr_of_biases,decomp_sigma=True)
  print(gauss_results)

def risingEdgeDiscretisation(files,trees,ch,total_number_channels):
  print("Rising edge discretisation analysing temporal distribution between CFD @ 20% and 50%")
  vars = []
  cut_conds = []
  toaThreshold = 4
  # 1 = 20%, 4 = 50%
  for ch_it in range(total_number_channels):
    var = "cfd["+str(ch_it)+"][1]-cfd["+str(ch_it)+"][4]"
    cut_cond = "event>-1"
    vars.append(var)
    cut_conds.append(cut_cond)

  nBins = 1000
  xLower = -6.0
  xUpper = 2.0

  arr_of_hists = []
  arr_of_biases = []

  for i in range(len(trees)):
    for j in range(len(vars)):
      bias = getBias(files[i])
      thisHist = hist_tree_file_basics(trees[i],files[i],vars[j],j,nBins,xLower,xUpper,bias,cut_conds[j],j,toaThreshold)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of time resolution (CFD@20%-CFD@50%) distribution", 800, 600)
  c1.SetLogy()
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(1, max_y)
  arr_of_hists[0].SetTitle("Comparison of time resolution (CFD@20%-CFD@50%) distribution")
  arr_of_hists[0].Draw("")
  for hist_to_draw in arr_of_hists[1:]:
    hist_to_draw.Draw("SAME")

  '''
  arr_of_fits = []
  for i in range(len(arr_of_hists)):
    thisFit = plot_fit_curves(xLower,xUpper,"gaus",arr_of_hists[i],i,arr_of_biases[i])
    thisFit.Draw("SAME")
    arr_of_fits.append(thisFit)
  '''

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists)):
    if ch == -1: ch_num = i+1
    if ch != -1: ch_num = ch+1
    legend.AddEntry(arr_of_hists[i], arr_of_biases[i] + " CH " + str(ch_num), "l")

  legend.Draw()
  c1.SaveAs("rising_edge_disc.png")


def waveform(files,trees,ch,total_number_channels):
  print("Waveform analysis")
  t_data = [[] for _ in range(total_number_channels)]
  w_data = [[] for _ in range(total_number_channels)]
  event_limit = 24500
  event_start = 4500
  for i in range(len(trees)):
    tree = trees[i]
    if ch != -1:
      event_cycle = 0
      for entry_index in range(event_start, event_limit+event_start):
        tree.GetEntry(entry_index)
        if (event_cycle > event_limit): continue
        event_cycle += 1
        t = tree.t
        w = tree.w
        pmax = tree.pmax
        if 0 < pmax[ch] < 200:
          t_data[0].extend(t[ch])
          w_data[0].extend(w[ch])
    else:
      for ch_it in range(total_number_channels):
        print("Channel number:")
        print(ch_it)
        event_cycle = 0
        for entry_index in range(event_start, event_limit+event_start):
          tree.GetEntry(entry_index)
          if (event_cycle > event_limit): continue
          event_cycle += 1
          t = tree.t
          w = tree.w
          pmax = tree.pmax
          if 0 < pmax[ch_it] < 200:
            t_data[ch_it].extend(t[ch_it])
            w_data[ch_it].extend(w[ch_it])
    
  t_data = [np.array(channel) for channel in t_data]
  w_data = [np.array(channel) for channel in w_data]

  fig, ax = plt.subplots(figsize=(10, 6))

  colors = plt.cm.viridis(np.linspace(0, 1, total_number_channels))

  for i in range(3):
    ax.scatter(t_data[i], w_data[i], label=f'Channel {i+1}', color=colors[i], s=1)

  ax.set_xlim(-0.5e-8,1e-8)
  ax.set_ylim(min(min(inner_array) for inner_array in w_data),max(max(inner_array) for inner_array in w_data))
  ax.set_xlabel('t / 10 ns', fontsize=14)
  ax.set_ylabel('Amplitude / V', fontsize=14)
  ax.set_title(f'Waveform for events between Ev No.{event_start} and Ev No.{event_limit}', fontsize=16)
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
      cut_cond = "event>-1 && negpmax["+str(j)+"] > -15 && pmax["+str(j)+"] > 6.0 && pmax["+str(j)+"] < 200"
      thisHist = hist_tree_file_basics(trees[i],files[i],"charge",j,nBins,xLower,xUpper,bias,cut_cond,j,0)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  arr_of_fits = []
  arr_of_divided_hists = []
  for i, hist in enumerate(arr_of_hists):
    hist_divided = hist.Clone(f"hist_divided_{i}")
    if i == 0:
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
      cut_cond = "event>-1 && negpmax["+str(j)+"] > -15 && pmax[" + str(j) + "] > 6.0 && pmax[" + str(j) + "] < 200"
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
    cut_cond = "event>-1 && negpmax["+str(ch_it)+"] > -15 && negpmax["+str(ch_it)+"] < 5 && pmax["+str(ch_it)+"] > 6.0 && pmax["+str(ch_it)+"] < 120 && negpmax["+str(ch_it+1)+"] > -15 && negpmax["+str(ch_it+1)+"] < 5 && pmax["+str(ch_it+1)+"] > 6.0 && pmax["+str(ch_it+1)+"] < 120 && cfd[0]["+str(toaThreshold)+"]>-2.2 && cfd[1]["+str(toaThreshold)+"]>-2.5 && cfd[2]["+str(toaThreshold)+"]>-2.0"
    vars.append(var)
    cut_conds.append(cut_cond)
  vars.append("cfd[0]["+str(toaThreshold)+"]-cfd["+str(total_number_channels-1)+"]["+str(toaThreshold)+"]")
  cut_conds.append("event>-1 && negpmax["+str(total_number_channels-1)+"] > -15 && negpmax["+str(total_number_channels-1)+"] < 5 && pmax["+str(total_number_channels-1)+"] > 6.0 && pmax["+str(total_number_channels-1)+"] < 120 && negpmax[0] > -15 && negpmax[0] < 5 && pmax[0] > 6.0 && pmax[0] < 120 && cfd[0]["+str(toaThreshold)+"]>-2.2 && cfd[1]["+str(toaThreshold)+"]>-2.5 && cfd[2]["+str(toaThreshold)+"]>-2.0")
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
  parser.add_argument('--doDiscretisation', action='store_true', help='Enable the doDiscretisation of the CFD rising edge option')
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
  if args.doDiscretisation: risingEdgeDiscretisation(file_array,tree_array,args.ch-1,total_number_channels)
  if args.doWaveform: waveform(file_array,tree_array,args.ch-1,total_number_channels)
  if args.csvOut: makeCSV(file_array,tree_array,total_number_channels)

if __name__ == "__main__":
    main()
