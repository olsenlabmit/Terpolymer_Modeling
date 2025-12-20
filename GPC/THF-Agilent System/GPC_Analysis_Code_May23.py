# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 15:14:32 2022

@author: ChemeGrad2020
"""

import logging

import numpy as np
import pandas as pd

import chem_analysis as chem
import chem_analysis.algorithms.baseline_correction as chem_bc

chem.logger_analysis.setLevel(logging.CRITICAL)


def main():
    # load data
    name='TOA50_TPD25_OTD25'
    file_name = r"C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/GPC/New_exports_24/"+str(name)+".csv"
    df = pd.read_csv(file_name, header=0, index_col=0)
    df = df.iloc[14:1600, :10]
    df.columns = ["LS1", "LS2", "LS3", "LS4", "LS5", "LS6", "LS7", "LS8", "UV", "RI"]
    print(df)
    # df=df.iloc[:6170,:2]
    # df.columns=["UV"]
    df.index.names = ["time (min)"]
    df = df.apply(pd.to_numeric, errors='coerce')
    print(df["RI"])
    # define calibration
    def cal_func(time: np.ndarray):
        return 10 ** (  - 0.4569*time + 10.663)
    
    # def cal_func_UV(time: np.ndarray):
    #     return 10 ** (0.00730 * (time**2) - 0.67753*time + 12.36205)
    
    # def cal_func_LS7(time: np.ndarray):#not updated april 2023
    #     return 10 ** (-0.467231*time + 10.932043)
        

    cal = chem.Cal(cal_func, lb=160, ub=1_090_000, name="calibration")
    # cal_UV = chem.Cal(cal_func_UV, lb=160, ub=1_090_000, name="UV calibration")
    # cal_LS7 = chem.Cal(cal_func_LS7, lb=160, ub=1_090_00, name= "LS7 calibration")

    # define signal
    signal = chem.SECSignal(ser=df["RI"], cal=cal, x_label="retention time (min)", y_label="signal")

    # data processing
    signal.pipeline.add(chem_bc.adaptive_polynomial_baseline)  # baseline correction
    signal.peak_picking(lb=10.8, ub=17)  # !!!!!!!!!!!! lb, and ub are lower and upper integration bounds for your peak (in retention time minutes) !!!!!!!!!!!!!!!!!!!!!!
    # print stats / plotting
    signal.stats(num_sig_figs=4)
    signal.plot()
   # signal.save("C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/GPC/"+str(name)+".png")

    #print("hi")


if __name__ == '__main__':
    main()