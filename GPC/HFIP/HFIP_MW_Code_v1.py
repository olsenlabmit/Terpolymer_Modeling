import logging
import os
import numpy as np
import pandas as pd
import chem_analysis as chem
import chem_analysis.algorithms.baseline_correction as chem_bc
import matplotlib.pyplot as plt
from scipy import stats

# Set logging level to CRITICAL
chem.logger_analysis.setLevel(logging.CRITICAL)

def process_file(file_path, output_dir):
    # Load data
    df = pd.read_csv(file_path, header=1, index_col=0)
    df = df.iloc[:30001, :2]
    df.columns = ["Y:"]
    df.index.names = ["X:"]
    df = df.apply(pd.to_numeric, errors='coerce')

    # Define calibration
    time = [6.922,8.729,10.623,
            6.711,8.067,9.998,
            7.490,9.311,10.992,
            ]
    MP_red_ri = [1_226_000, 130_000, 13_630, 
                 1_795_000, 273_600, 33_240,
                 536_500, 74000,6960
                 ]
        # Log-transform the refractive index data (y-data)
    log_MPR = np.log(MP_red_ri)
    
    # Perform linear regression on the log-transformed data (log(RI) vs time)
    slope, intercept, r_value, p_value, std_err = stats.linregress(time, log_MPR)
    
    # Create the linear calibration function
    def cal_func_RI(t):
        return np.exp(slope * t + intercept)
    
    # Generate a range of time values for plotting the linear fit
    time_range = np.linspace(min(time)-0.5, max(time)+0.5, 5000)
    plt.scatter(time, log_MPR, color='red', label='Original Data', zorder=5)
    plt.plot(time_range, slope * time_range + intercept, color='blue', label='Linear Fit')
    plt.title("Linear Fit for Refractive Index Calibration")
    plt.xlabel("Time")
    plt.ylabel("Log of Refractive Index (log(RI))")
    plt.legend()
    plt.show()

    cal_RI = chem.Cal(cal_func_RI, lb=200, ub=2_000_000, name="RI calibration")
    signal = chem.SECSignal(ser=df["Y:"], cal=cal_RI, x_label="X:", y_label="signal")

    # Data processing
    signal.pipeline.add(chem_bc.adaptive_polynomial_baseline)  # baseline correction
    signal.peak_picking(lb=6.8, ub=12.2)  # Set the lower and upper bounds for peak picking

    # Create output subdirectory for this file
    file_name = os.path.basename(file_path)
    output_subdir = os.path.join(output_dir, os.path.splitext(file_name)[0])
    os.makedirs(output_subdir, exist_ok=True)

    # Save stats to CSV
    stats_file = os.path.join(output_subdir, f"{os.path.splitext(file_name)[0]}_stats.csv")
    signal.stats(num_sig_figs=4)
    # signal.stats_df.to_csv(stats_file)

    # Save plot to file
    plot_file = os.path.join(output_subdir, f"{os.path.splitext(file_name)[0]}_plot.png")
    signal.plot()
    plt.savefig(plot_file)
    plt.close()

def process_files_in_range(input_dir, output_dir, start_sample, end_sample):
    # Get list of all CSV files in the input directory
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    
    # Sort files by name (optional, depending on how the files are named)
    csv_files.sort()

    # Filter files based on the sample range (e.g., start at sample 1, go for sample 3)
    selected_files = csv_files[start_sample - 1:end_sample]
    
    # Print the name of the folder (directory) being processed
    folder_name = os.path.basename(input_dir)
    print(f"Processing files in folder: {folder_name}")
    
    # Process each file in the range
    for csv_file in selected_files:
        file_path = os.path.join(input_dir, csv_file)
        print(f"Processing file: {csv_file}___________________________")  # Optionally print each file being processed
        process_file(file_path, output_dir)

def main():
    # Input and output directories
    #input_dir = r"C:\Users\nmamr\OneDrive - Massachusetts Institute of Technology\Graduate School\Polymer Science Research\Furanoate Library Paper\Molecular Weight\data"
    #output_dir = r"C:\Users\nmamr\OneDrive - Massachusetts Institute of Technology\Graduate School\Polymer Science Research\Furanoate Library Paper\Molecular Weight\analysis"
    input_dir = r"C:\Users\ChemeGrad2020\Documents\Grad School\Research\Ter-Polymer_Library_Biodeg\GPC\HFIP-Revisions"
    output_dir = r"C:\Users\ChemeGrad2020\Documents\Grad School\Research\Ter-Polymer_Library_Biodeg\GPC\HFIP-Revisions\analysis"

    # Define the sample range (e.g., start from sample 1 to sample 3)
    start_sample = 34
    end_sample = start_sample

    # Process the files in the specified range
    process_files_in_range(input_dir, output_dir, start_sample, end_sample)

if __name__ == '__main__':
    main()
