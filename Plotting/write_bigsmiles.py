# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 15:59:07 2024

@author: ChemeGrad2020
"""

#write bigsmiles
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
import csv
import random
import pandas as pd
import math

terpolymers = pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Library_Data_Final.csv', header=0)
julia=pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Library_Data_Final_Julia.csv', header=0,sep='\t')
monomers = pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Monomer Info_DMTDMI.csv')

mon_names=monomers['Polymer Name Abbr.'].tolist()
smilecont=monomers['SMILES for Bigsmiles'].tolist()
names=julia['Polymer_Name_P'].tolist()
print(terpolymers['Polymer_Name_P'])
print(julia['Polymer_Name_P'])
with open('namesandbigsmiles_julia.csv', mode='w') as file:
    file_writer=csv.writer(file, delimiter=',')
    file_writer.writerow(['Polymer_Name_P','BIGSMILES'])

for i in range(len(names)):
    string='{[]'
    string=string + smilecont[mon_names.index(names[i][0:3])] +','
    string=string + smilecont[mon_names.index(names[i][6:9])] +','
    string=string + smilecont[mon_names.index(names[i][12:15])] +'[]}'
    with open('namesandbigsmiles_julia.csv', mode='a') as file:
        file_writer=csv.writer(file, delimiter=',')
        file_writer.writerow([names[i],string])
        
    