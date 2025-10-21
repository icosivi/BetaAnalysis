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
import math

#dir_name = "../TB_SPS_June_CMS_sensors/Run0_TB_SPS_LFoundry-K1/"

var_dict = {"tmax":"t_{max} / ns" , "pmax":"p_max / mV" , "negpmax":"-p_max / mV", "charge":"Q / fC", "area_new":"Area / pWb" , "rms":"RMS / mV"}

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

def get_fit_results(arr_of_results_to_fit,arr_of_biases,decomp_sigma=False,simplified=False):
  arr_of_mean = []
  arr_of_sigma = []
  arr_of_ampl = []
  #arr_of_chi2 = []
  #arr_of_ndf = []
  arr_of_red_chi2 = []
  #arr_of_prob = []

  for fit_func in arr_of_results_to_fit:
    mean = fit_func.GetParameter(1)  # Mean of the gauss distribution
    sigma = fit_func.GetParameter(2) # Sigma of the gauss distribution
    amplitude = fit_func.GetParameter(0)  # Amplitude of the gauss distribution
    chi2 = fit_func.GetChisquare()  # Chi-squared value of the fit
    ndf = fit_func.GetNDF()  # Number of degrees of freedom
    #prob = fit_func.GetProb()  # Probability of the fit result
    arr_of_mean.append(round_to_sig_figs(mean,4))
    arr_of_sigma.append(round_to_sig_figs(sigma,4))
    arr_of_ampl.append(round_to_sig_figs(amplitude,3))
    #arr_of_chi2.append(chi2)
    #arr_of_ndf.append(ndf)
    print(ndf)
    if ndf == 0:
      ndf = 1
    arr_of_red_chi2.append(round_to_sig_figs((chi2/ndf),3))
    #arr_of_prob.append(prob)

  df_of_results = pd.DataFrame({
    "Bias": arr_of_biases,
    "Mean": arr_of_mean,
    "Sigma": arr_of_sigma,
    "Amplitude": arr_of_ampl,
    #"Chi2": arr_of_chi2,
    #"NDF": arr_of_ndf,
    "RChi2": arr_of_red_chi2
    #"Prob": arr_of_prob
  })

  if decomp_sigma:
    if len(arr_of_sigma) < 3:
      sig1 = np.sqrt(arr_of_sigma[0]**2 - 0.015**2)
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
  elif var == "charge":
    thisHist = root.TH1F("hist"+biasVal, "area_new;"+var_dict[var]+";Events", nBins, xLower, xUpper)
    tree.Draw("area_new["+str(ch)+"]/4.7>>hist"+biasVal,cut_cond)
  elif var == "width":
    thisHist = root.TH1F("hist"+biasVal, "width[0][1];Width / ns;Events", nBins, xLower, xUpper)
    tree.Draw("width[0]["+str(toaThreshold)+"]>>hist"+biasVal,cut_cond)
  else:
    thisHist = root.TH1F("hist"+biasVal, var+";"+var_dict[var]+";Events", nBins, xLower, xUpper)
    tree.Draw(var+"["+str(ch)+"]>>hist"+biasVal,cut_cond)
  thisHist.SetLineWidth(2)
  thisHist.SetLineColor(index+1+ch)
  return thisHist

def plot_fit_curves(xLower,xUpper,fit_type,hist_to_fit,index,biasVal):
  thisFit = TF1(fit_type+"_hist"+biasVal, fit_type, xLower, xUpper)
  hist_to_fit.Fit(thisFit, "Q")
  thisFit.SetLineWidth(3)
  thisFit.SetLineColor(index+1)
  #thisFit.SetLineStyle(2)
  return thisFit

def waveform(files,trees,ch,total_number_channels):
  t_data = [[] for _ in range(total_number_channels)]
  w_data = [[] for _ in range(total_number_channels)]
  event_limit = 51000

  for i in range(len(trees)):
    tree = trees[i]
    event_cycle = 50000
    for entry in tree:
      if (event_cycle > event_limit): continue
      event_cycle += 1
      t = entry.t
      w = entry.w
      pmax = entry.pmax
      if 0 < pmax[0] < 220:
        t_data[i].extend(t[0])
        w_data[i].extend(w[0])

  t_data = [np.array(channel) for channel in t_data]
  w_data = [np.array(channel) for channel in w_data]

  fig, ax = plt.subplots(figsize=(10, 6))

  colors = plt.cm.viridis(np.linspace(0, 1, total_number_channels))

  ax.scatter(t_data[0], w_data[0], label=f'W5 Preirrad', color=colors[0])
  ax.scatter(t_data[1], w_data[1], label=f'W13 Irrad', color=colors[1])

  ax.set_xlim(0.8e-8,1.4e-8)
  ax.set_ylim(min(min(inner_array) for inner_array in w_data),max(max(inner_array) for inner_array in w_data))
  ax.set_xlabel('t / ns', fontsize=14)
  ax.set_ylabel('Amplitude / V', fontsize=14)
  ax.set_title(f'Waveform comparison, first {event_limit} events per channel', fontsize=16)
  ax.legend(loc = "upper right")
  ax.grid(True)
  plt.tight_layout()
  plt.savefig("waveform_comparison.png",facecolor='w')

  var = "width"
  nBins_w = 200
  arr_of_hists = []
  arr_of_biases = []
  cut_cond = "event>-1 && negpmax[0] > -15 && negpmax[0] < 5 && pmax[1] < 120 && pmax[0] > 10.0 && negpmax[1] > -18 && pmax[0] < 120.0"

  xLower_w = 0
  xUpper_w = 40

  for i in range(len(trees)):
    bias = getBias(files[i])
    thisHist = hist_tree_file_basics(trees[i],files[i],var,0,nBins_w,xLower_w,xUpper_w,bias,cut_cond,i,1)
    arr_of_hists.append(thisHist)
    arr_of_biases.append(bias)

  c1 = root.TCanvas("c1", "Comparison of width distribution", 800, 600)
  c1.SetLogy()
  max_y_w = max(hist.GetMaximum() for hist in arr_of_hists) * 1.05
  arr_of_hists[0].GetYaxis().SetRangeUser(1, max_y_w)
  arr_of_hists[0].SetTitle("Comparison of width distribution")
  arr_of_hists[0].Draw()
  for hist_to_draw in arr_of_hists:
    hist_to_draw.Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  wafer = ["W5 Preirrad","W13 35e14"]
  for i in range(len(arr_of_hists)):
    legend.AddEntry(arr_of_hists[i], wafer[i], "l")

  legend.Draw()

  c1.SaveAs("comp_width_test.png")

  var_a = "area_new"
  nBins_a = 500
  arr_of_hists2 = []
  arr_of_biases2 = []

  xLower_a = 0
  xUpper_a = 100

  for i in range(len(trees)):
    bias = getBias(files[i])
    thisHist = hist_tree_file_basics(trees[i],files[i],var_a,0,nBins_a,xLower_a,xUpper_a,bias,cut_cond,i,1)
    arr_of_hists2.append(thisHist)
    arr_of_biases2.append(bias)

  arr_of_divided_hists2 = []
  for i, hist in enumerate(arr_of_hists2):
    hist_divided2 = hist.Clone(f"hist_divided2_{i}")
    for bin in range(1, hist_divided2.GetNbinsX() + 1):
      bin_content = hist_divided2.GetBinContent(bin)
      #hist_divided2.SetBinContent(bin, bin_content / 5.0) # Mignone + 40 dB
      hist_divided2.SetBinContent(bin, bin_content / 4.7) # SC + 20 dB
    arr_of_divided_hists2.append(hist_divided2)

  c2 = root.TCanvas("c1", "Comparison of area distribution", 800, 600)
  c2.SetLogy()
  max_y_a = max(hist.GetMaximum() for hist in arr_of_hists2) * 1.05
  arr_of_hists2[0].GetYaxis().SetRangeUser(1, max_y_a)
  arr_of_hists2[0].SetTitle("Comparison of area distribution")
  arr_of_hists2[0].Draw()
  for hist_to_draw in arr_of_divided_hists2:
    hist_to_draw.Draw("SAME")

  legend = root.TLegend(0.7, 0.7, 0.9, 0.9)
  for i in range(len(arr_of_hists2)):
    legend.AddEntry(arr_of_hists2[i], wafer[i], "l")

  legend.Draw()

  c2.SaveAs("comp_charge_test.png")

  # binsw 200, binsa 500
  hist2D1 = root.TH2F("hist2D1", wafer[0], 30, 0.5, 2.0, 32, 0, 16)
  
  for entry in trees[0]:
    width_value = entry.width[0][1] if len(entry.width) > 0 and len(entry.width[0]) > 1 else np.nan
    area_value = entry.area_new[0] / 4.7 if len(entry.area_new) > 0 else np.nan # SC + 20 dB
    #area_value = entry.area_new[0] / 5 if len(entry.area_new) > 0 else np.nan # Mignone + 40 dB
    if entry.pmax[0] > 10.0 and entry.pmax[0] < 120.0 and entry.pmax[1] < 120.0 and entry.cfd[1][1]<2 and entry.cfd[1][1]>-2 and entry.cfd[0][1]>10:
      hist2D1.Fill(width_value, area_value)

  hist2D1.Draw("COLZ")
  hist2D1.GetXaxis().SetTitle("Width / ns")
  hist2D1.GetYaxis().SetTitle("Charge / fC")
  hist2D1.SetTitle(wafer[0])

  c2D1 = root.TCanvas("c2D1", "W5", 800, 600)
  hist2D1.Draw("COLZ")
  c2D1.SaveAs("area_v_width_W5.png")

  hist2D2 = root.TH2F("hist2D2", wafer[1], 30, 0.6, 2.1, 32, 0, 14)
  
  for entry in trees[1]:
    width_value = entry.width[0][1] if len(entry.width) > 0 and len(entry.width[0]) > 1 else np.nan
    area_value = entry.area_new[0] / 4.7 if len(entry.area_new) > 0 else np.nan # SC + 20 dB
    #area_value = entry.area_new[0] / 5 if len(entry.area_new) > 0 else np.nan # Mignone + 40 dB
    if entry.pmax[0] > 6.0 and entry.pmax[0] < 200.0 and entry.pmax[1] < 200.0 and entry.cfd[1][1]<2 and entry.cfd[1][1]>-2 and entry.cfd[0][1]>10:
      hist2D2.Fill(width_value, area_value)

  hist2D2.Draw("COLZ")
  hist2D2.GetXaxis().SetTitle("Width / ns")
  hist2D2.GetYaxis().SetTitle("Charge / fC")
  hist2D2.SetTitle(wafer[1])

  c2D2 = root.TCanvas("c2D2", "W13", 800, 600)
  hist2D2.Draw("COLZ")
  c2D2.SaveAs("area_v_width_W13.png")

def main():
  parser = argparse.ArgumentParser(description='Read .root files into an array.')
  parser.add_argument('files', metavar='F', type=str, nargs='+',
                      help='List of .root files or wildcard pattern (*.root)')
  parser.add_argument('--ch', type=int, nargs='?', default=0,
                      help="For analysis of a specific channel, choose an integer in the range 1 to the total number of channels")
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


  if args.doWaveform: waveform(file_array,tree_array,args.ch-1,total_number_channels)

if __name__ == "__main__":
    main()
