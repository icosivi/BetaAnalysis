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

dir_name = "../TB_SPS_June_CMS_sensors/Run0_TB_SPS_LFoundry-K1/"

var_dict = {"tmax":"t_{max} / ns" , "pmax":"p_max / pWb" , "charge":"Q / fC", "area_new":"Area / pWb" , "rms":"RMS / mV"}

doTMAXTest = False
doPMAXTest = False
doLandauToArea = False
doRMSNoise = False
doTimeRes = True

'''
time res v bias
collected charge v bias
rms v bias

t1-t2 Cividec
t3-t4 Minicircuit

Analysis cuts pmax 30 to 300 V
'''

def get_fit_results(arr_of_results_to_fit,arr_of_biases):
  arr_of_mean = []
  arr_of_sigma = []
  arr_of_ampl = []
  arr_of_chi2 = []
  arr_of_ndf = []
  arr_of_red_chi2 = []
  arr_of_prob = []

  for fit_func in arr_of_results_to_fit:
    mean = fit_func.GetParameter(1)  # Mean of the gauss distribution
    sigma = fit_func.GetParameter(2) # Sigma of the gauss distribution
    amplitude = fit_func.GetParameter(0)  # Amplitude of the gauss distribution
    chi2 = fit_func.GetChisquare()  # Chi-squared value of the fit
    ndf = fit_func.GetNDF()  # Number of degrees of freedom
    prob = fit_func.GetProb()  # Probability of the fit result
    arr_of_mean.append(mean)
    arr_of_sigma.append(sigma)
    arr_of_ampl.append(amplitude)
    arr_of_chi2.append(chi2)
    arr_of_ndf.append(ndf)
    arr_of_red_chi2.append(chi2/ndf)
    arr_of_prob.append(prob)

  df_of_results = pd.DataFrame({
    "Bias": arr_of_biases,
    "Mean": arr_of_mean,
    "Sigma": arr_of_sigma,
    "Amplitude": arr_of_ampl,
    "Chi2": arr_of_chi2,
    "NDF": arr_of_ndf,
    "RChi2": arr_of_red_chi2,
    "Prob": arr_of_prob
  })

  return df_of_results

def getBias(filename):
  pattern = r"_(\d{2,3}V)."
  match = re.search(pattern, str(filename))
  if match:
    return str(match.group(1))
  else:
    return None

def hist_tree_file_basics(tree,file,var,index,nBins,xLower,xUpper,biasVal,cut_cond,ch):
  if var == "cfd["+str(ch)+"][1]-cfd["+str(ch+1)+"][1]":
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

def plot_fit_curves(xLower,xUpper,fit_type,hist_to_fit,index,bias):
  thisFit = TF1(fit_type+"_hist"+bias, fit_type, xLower, xUpper)
  hist_to_fit.Fit(thisFit, "Q")
  thisFit.SetLineColor(index+1)
  thisFit.SetLineWidth(3)
  return thisFit

def tMaxTest(files,trees,ch,total_number_channels):
  print("tmax test run")
  var = "tmax"
  nBins = 1000
  cut_cond = ""
  arr_of_hists = []
  arr_of_biases = []

  xLower = trees[0].GetMinimum(var)
  xUpper = trees[0].GetMaximum(var)

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      print(files[i])
      print(bias)
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j)
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

  xLower = 0
  xUpper = trees[0].GetMaximum(var)
  cut_cond = "pmax["+str(ch)+"] > 30 && pmax["+str(ch)+"] < 300"

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j)
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

def landauToArea(files,trees,ch,total_number_channels):
  print("Fitting Landau curve to area distribution")
  var = "area_new"
  nBins = 200
  xLower = 0
  xUpper = 300
  cut_cond = "pmax["+str(ch)+"] > 30 && pmax["+str(ch)+"] < 300"

  arr_of_hists = []
  arr_of_biases = []

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j)
        arr_of_hists.append(thisHist)
        arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of area_new distribution", 800, 600)
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(0, max_y)
  arr_of_hists[0].SetTitle("Comparison of area_new distribution")
  arr_of_hists[0].Draw()
  for hist_to_draw in arr_of_hists:
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

def landauToCharge(files, trees, ch,total_number_channels):
  print("Fitting Landau curve to charge distributions")
  var = "charge"
  nBins = 200
  xLower = 0
  xUpper = 60
  cut_cond = "pmax[" + str(ch) + "] > 30 && pmax[" + str(ch) + "] < 300"

  arr_of_hists = []
  arr_of_biases = []

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j)
        arr_of_hists.append(thisHist)
        arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of charge distributions", 800, 600)
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(0, max_y)
  arr_of_hists[0].SetTitle("Comparison of charge distribution")
  arr_of_hists[0].Draw()
  for hist_to_draw in arr_of_hists[1:]:
    hist_to_draw.Draw("SAME")

  arr_of_fits = []
  for i in range(len(arr_of_hists)):
    print("Fit results for this dist")
    hist_divided = arr_of_hists[i].Clone()

    for bin in range(1, hist_divided.GetNbinsX() + 1):
      x_val = hist_divided.GetXaxis().GetBinCenter(bin)
      hist_divided.SetBinContent(bin, hist_divided.GetBinContent(bin) / 4.7)

    thisFit = plot_fit_curves(xLower, xUpper, "landau", hist_divided, i, arr_of_biases[i])
    arr_of_fits.append(thisFit)
    thisFit.Draw("SAME")

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
  xUpper = trees[0].GetMaximum(var)
  cut_cond = "pmax["+str(ch)+"] > 30 && pmax["+str(ch)+"] < 300"

  arr_of_hists = []
  arr_of_biases = []

  if total_number_channels == 1:
    for i in range(len(trees)):
      bias = getBias(files[i])
      thisHist = hist_tree_file_basics(trees[i],files[i],var,i,nBins,xLower,xUpper,bias,cut_cond,ch)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)
  else:
    for i in range(len(trees)):
      for j in range(total_number_channels):
        bias = getBias(files[i])
        thisHist = hist_tree_file_basics(trees[i],files[i],var,j,nBins,xLower,xUpper,bias,cut_cond,j)
        arr_of_hists.append(thisHist)
        arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of RMS distribution", 800, 600)
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(0, max_y)
  arr_of_hists[0].SetTitle("Comparison of RMS distribution")
  arr_of_hists[0].Draw()
  for hist_to_draw in arr_of_hists:
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
  c1.SaveAs("Gauss_analysis.png")

  gauss_results = get_fit_results(arr_of_fits,arr_of_biases)
  print(gauss_results)

def timeRes(files,trees,ch,total_number_channels):
  print("Time resolution analysis")
  vars = []
  cut_conds = []
  if ch == -1:
    for ch_it in range(total_number_channels-1):
      var = "cfd["+str(ch_it)+"][1]-cfd["+str(ch_it+1)+"][1]"
      cut_cond = "pmax["+str(ch_it)+"] > 30 && pmax["+str(ch_it)+"] < 300"
      vars.append(var)
      cut_conds.append(cut_cond)
  else:
    var = "cfd["+str(ch)+"][1]-cfd["+str(ch+1)+"][1]"
    cut_cond = "pmax["+str(ch)+"] > 30 && pmax["+str(ch)+"] < 300"
    vars.append(var)
    cut_conds.append(cut_cond)

  nBins = 200
  xLower = 0
  xUpper = trees[0].GetMaximum(var)

  arr_of_hists = []
  arr_of_biases = []

  for i in range(len(trees)):
    for j in range(len(vars)):
      bias = getBias(files[i])
      thisHist = hist_tree_file_basics(trees[i],files[i],vars[j],j,nBins,xLower,xUpper,bias,cut_conds[j],j)
      arr_of_hists.append(thisHist)
      arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of time resolution (CFD@20%) distribution", 800, 600)
  max_y = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(0, max_y)
  arr_of_hists[0].SetTitle("Comparison of time resolution (CFD@20%) distribution")
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

  gauss_results = get_fit_results(arr_of_fits,arr_of_biases)
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
        if 30 < pmax[ch] < 300:
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
          if 30 < pmax[ch_it] < 300:
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

def main():
  parser = argparse.ArgumentParser(description='Read .root files into an array.')
  parser.add_argument('files', metavar='F', type=str, nargs='+',
                      help='List of .root files or wildcard pattern (*.root)')
  parser.add_argument('--ch', type=int, nargs='?', default=0,
                      help="For analysis of a specific channel, choose an integer in the range 1 to the total number of channels")
  parser.add_argument('--doTMAXTest', action='store_true', help='Enable the doTMAXTest option')
  parser.add_argument('--doPMAXTest', action='store_true', help='Enable the doPMAXTest option')
  parser.add_argument('--doLandauToArea', action='store_true', help='Enable the doLandauToArea option')
  parser.add_argument('--doChargeDist', action='store_true', help='Enable the doChargeDist option')
  parser.add_argument('--doRMSNoise', action='store_true', help='Enable the doRMSNoise option')
  parser.add_argument('--doTimeRes', action='store_true', help='Enable the doTimeRes option')
  parser.add_argument('--doWaveform', action='store_true', help='Enable the doWaveform option')
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
  if args.doLandauToArea: landauToArea(file_array,tree_array,args.ch-1,total_number_channels)
  if args.doChargeDist: landauToCharge(file_array,tree_array,args.ch-1,total_number_channels)
  if args.doRMSNoise: rmsNoise(file_array,tree_array,args.ch-1,total_number_channels)
  if args.doTimeRes: timeRes(file_array,tree_array,args.ch-1,total_number_channels)
  if args.doWaveform: waveform(file_array,tree_array,args.ch-1,total_number_channels)

if __name__ == "__main__":
    main()
