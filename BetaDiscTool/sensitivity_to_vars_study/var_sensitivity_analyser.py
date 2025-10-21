import ROOT
import numpy as np
import math
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt


file = ROOT.TFile.Open("stats_Sr_Run5_Ch1-180V_Ch2-130V_Ch3-2750V_trig230V.root")
tree = file.Get("Analysis")

cfd_indices = list(range(5))
width_indices = list(range(5))

vars_to_analyze = [
  'pmax', 'negpmax', 'tmax', 'negtmax', 'area_new',
  'uarea_new', 'risetime', 'falltime', 'dvdt',
  'dvdt_2080', 'tot', 'rms'
]

for i in width_indices:
  vars_to_analyze.append(f'width_{(i+1)*10}pc')

pmax_cuts = {1: 10, 2: 12}

data_ch = {
    1: {var: [] for var in vars_to_analyze},
    2: {var: [] for var in vars_to_analyze}
}
labels_ch = {1: [], 2: []}
n_entries = tree.GetEntries()

for entry in range(n_entries):
  tree.GetEntry(entry)
  pmax = list(tree.pmax)

  for ch in [1, 2]:
    if ch >= len(pmax):
      continue

    label = 1 if pmax[ch] > pmax_cuts[ch] else 0
    labels_ch[ch].append(label)

    for var in [
      'pmax', 'negpmax', 'tmax', 'negtmax',
      'area_new', 'uarea_new', 'risetime',
      'falltime', 'dvdt', 'dvdt_2080', 'tot', 'rms'
    ]:
      val = getattr(tree, var)
      data_ch[ch][var].append(val[ch] if ch < len(val) else np.nan)

    cfd_vec = tree.cfd[ch] if ch < len(tree.cfd) else []
    width_vec = tree.width[ch] if ch < len(tree.width) else []

    for i in width_indices:
      varname = f'width_{(i+1)*10}pc'
      val = width_vec[i] if i < len(width_vec) else np.nan
      data_ch[ch][varname].append(val)

roc_results = {1: {}, 2: {}}

print("\nVariable Sensitivities (ROC AUC):")
for var in vars_to_analyze:
  for ch in [1, 2]:
    values = np.array(data_ch[ch][var])
    ch_labels = np.array(labels_ch[ch])

    mask = (
      ~np.isnan(values) &
      ~np.isinf(values) &
      (np.abs(values) < 1e10)
    )
    clean_labels = ch_labels[mask]
    clean_values = values[mask]

    if len(np.unique(clean_labels)) < 2 or len(clean_values) < 10:
      print(f"Ch{ch} — {var:20s} : Skipped (insufficient data)")
      continue

    try:
      auc = roc_auc_score(clean_labels, clean_values)
      roc_results[ch][var] = auc
      print(f"Ch{ch} — {var:20s} : {auc:.3f}")
    except Exception as e:
      print(f"Ch{ch} — {var:20s} : Error calculating AUC → {str(e)}")

rename_dict = {
  "pmax": r"\mathrm{Ampl.}^{\mathrm{+ve}}",
  "negpmax": r"\mathrm{Ampl.}^{\mathrm{-ve}}",
  "tmax": r"t_{\mathrm{ampl.^{+ve}}}",
  "negtmax": r"t_{\mathrm{ampl.^{-ve}}}",
  "area_new": r"\mathrm{Area^{+ve}}",
  "uarea_new": r"\mathrm{Area^{-ve}}",
  "risetime": r"t_{\mathrm{rise}}",
  "falltime": r"t_{\mathrm{fall}}",
  "dvdt": r"\mathrm{Slew~rate}~\left. \frac{dV}{dt} \right|_{20\%}",
  "dvdt_2080": r"\mathrm{Slew~rate}~\left. \frac{dV}{dt} \right|^{80\%}_{20\%}",
  "tot": r"t_{\mathrm{total}}",
  "rms": r"\mathrm{RMS~noise}",
  "width_10pc": r"\mathrm{Width_{10\%~ampl.}}",
  "width_20pc": r"\mathrm{Width_{20\%~ampl.}}",
  "width_30pc": r"\mathrm{Width_{30\%~ampl.}}",
  "width_40pc": r"\mathrm{Width_{40\%~ampl.}}",
  "width_50pc": r"\mathrm{Width_{50\%~ampl.}}"
}

roc_results_renamed = {
  1: {rename_dict.get(k, k) or k: v for k, v in roc_results[1].items()},
  2: {rename_dict.get(k, k) or k: v for k, v in roc_results[2].items()}
}

all_vars = list(rename_dict.get(k, k) or k for k in vars_to_analyze)

auc1 = [roc_results_renamed[1].get(v, np.nan) for v in all_vars]
auc2 = [roc_results_renamed[2].get(v, np.nan) for v in all_vars]

latex_labels = [f"${v}$" for v in all_vars]
y = np.arange(len(all_vars))
bar_width = 0.4

plt.figure(figsize=(10, max(6, len(latex_labels) * 0.35)))
plt.barh(y - bar_width/2, auc1, height=bar_width, color='dodgerblue', label='20 μm (Ch 1) 180 V')
plt.barh(y + bar_width/2, auc2, height=bar_width, color='red', label='20 μm (Ch 2) 130 V')

plt.yticks(y, latex_labels)
plt.xlabel("ROC AUC")
plt.legend()
plt.gca().invert_yaxis()
plt.tight_layout()
plt.grid(True, axis='x', linestyle='--', alpha=0.6)
plt.savefig("sensitivity_plot.png", dpi=300, bbox_inches='tight', facecolor='white')
