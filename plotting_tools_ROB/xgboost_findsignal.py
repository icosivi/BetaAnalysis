import ROOT
import pandas as pd
import numpy as np
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
import argparse

# Read branches from ROOT files to train on
def read_root_files(file_list):
    data_frames = []
    
    for file in file_list:
        f = ROOT.TFile.Open(file)
        tree = f.Get("Analysis")
        
        pmax = []
        rms = []
        charge = []
        w = []
        t = []
        negpmax = []
        cfd = []

        for event in tree:
            pmax.append(event.pmax[0])
            rms.append(event.rms[0])
            charge.append(event.area_new[0])
            #w.append(event.w[0])
            #t.append(event.t[0])
            negpmax.append(event.negpmax[0])
            cfd.append(event.cfd[0][1])
            #current.append(event.I[0])
            #voltage.append(event.V[0])

        df = pd.DataFrame({
            'pmax': pmax,
            'rms': rms,
            'charge': charge,
            'negpmax': negpmax,
            'cfd': cfd
        })
        data_frames.append(df)

    return pd.concat(data_frames, ignore_index=True)

def prepare_training_data(df):
    # cut on pmax and label signal (1) and background (0)
    pmax_cut = 35
    df['label'] = np.where(df['pmax'] > pmax_cut, 1, 0)
    
    X = df.drop(columns=['pmax', 'label'])
    y = df['label']
    
    return X, y

def train_classifier(X, y):
    # Split data into training and validation sets
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.3, random_state=42)
    
    # Initialize XGBoost classifier
    clf = XGBClassifier(use_label_encoder=False, eval_metric="logloss", learning_rate=0.005, n_estimators=5)
    '''
    clf = XGBClassifier(
        use_label_encoder=False,
        eval_metric="logloss",
        n_estimators=5,
        learning_rate=0.005,
        max_depth=6,
        subsample=0.8,
        colsample_bytree=0.8
    )
    '''
    
    # Train the model
    clf.fit(X_train, y_train)
    
    # Evaluate the model using ROC AUC score
    y_val_pred = clf.predict_proba(X_val)[:, 1]
    auc_score = roc_auc_score(y_val, y_val_pred)
    print(f"Validation AUC: {auc_score}")
    
    return clf

# Function to evaluate the model and predict signal probability on new data
def evaluate_model(clf, eval_file):
    # Read evaluation data
    df_eval = read_root_files([eval_file])
    
    # Use all features except pmax for evaluation
    X_eval = df_eval.drop(columns=["pmax"])
    
    # Predict signal probability
    df_eval['signal_prob'] = clf.predict_proba(X_eval)[:, 1]
    
    # Print or save results
    print(df_eval[['pmax', 'rms', 'charge', 'negpmax', 'cfd', 'signal_prob']])
    
    # Save the result to a CSV file
    output_file = "evaluation_results.csv"
    df_eval.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Signal and background classifier using specific branches from ROOT files")
    parser.add_argument("--train_files", nargs='+', help="List of ROOT files for training", required=True)
    parser.add_argument("--eval_file", help="ROOT file to evaluate the model", required=True)
    
    args = parser.parse_args()
    
    train_files = args.train_files
    print(f"[FINDSIGNAL] : Reading training files: {train_files}")
    df_train = read_root_files(train_files)
    
    X_train, y_train = prepare_training_data(df_train)
    
    print("[FINDSIGNAL] : Training XGBoost algorithm with following parameters: ")
    clf = train_classifier(X_train, y_train)
    
    eval_file = args.eval_file
    print(f"[FINDSIGNAL] : Evaluating the model on: {eval_file}")
    evaluate_model(clf, eval_file)

if __name__ == "__main__":
    main()
