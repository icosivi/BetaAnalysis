import ROOT
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy.stats import ks_2samp


# Function to read real pmax distribution from eval_file.root
def read_real_pmax(file):
    f = ROOT.TFile.Open(file)
    tree = f.Get("Analysis")  # Ensure this is the correct TTree name
    
    # Prepare list to hold pmax values
    pmax_real = []
    
    # Loop over the entries in the tree and collect pmax values
    for event in tree:
        pmax_real.append(event.pmax[0])  # Assuming pmax is a scalar
    
    f.Close()
    
    return np.array(pmax_real)

# Function to compare real and reconstructed pmax distributions
def compare_distributions(eval_file, results_file):
    # Step 1: Read the real pmax distribution from the ROOT file
    pmax_real = read_real_pmax(eval_file)
    
    # Step 2: Read the evaluation results (reconstructed data)
    df_eval = pd.read_csv(results_file)
    
    # Step 3: Reconstruct pmax based on signal probability (weight the distribution by signal_prob)
    pmax_reconstructed = df_eval['pmax']
    signal_prob = df_eval['signal_prob']
    
    # You can use signal_prob as a weight for the reconstructed pmax values.
    weights = signal_prob
    
    # Step 4: Plot the real and reconstructed pmax distributions
    plt.figure(figsize=(10, 6))
    
    # Real pmax distribution
    plt.hist(pmax_real, bins=50, alpha=0.5, label='Real pmax', density=True, color='blue')
    
    # Labels and legend
    plt.xlabel('pmax')
    plt.ylabel('Density')
    plt.title('Real vs. Reconstructed pmax Distribution')

    ks_stat, ks_p_value = ks_2samp(pmax_real, pmax_reconstructed)
    # Reconstructed pmax distribution (weighted by signal_prob)
    plt.hist(pmax_reconstructed, bins=50, alpha=0.5, label=f'Reconstructed pmax (weighted)\nS_KS = {ks_stat}, p_KS = {ks_p_value}', 
             weights=weights, density=True, color='red')

    plt.legend()
    
    # Show the plot
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Compare real and reconstructed pmax distributions")
    
    # Arguments for input files
    parser.add_argument("--eval_file", help="Original ROOT file to compare to (e.g., eval_file.root)", required=True)
    parser.add_argument("--results_file", help="Evaluation results CSV file (e.g., evaluation_results.csv)", required=True)
    
    args = parser.parse_args()
    
    # Compare distributions using the provided files
    compare_distributions(args.eval_file, args.results_file)

if __name__ == "__main__":
    main()
