# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 13:47:15 2024

@author: ChemeGrad2020
"""

#make abbreviation names for each of the original library polymers to allow for string searching
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
import csv
import random
import pandas as pd



polymers = pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Original_Library_Degradability.csv', header=0)
print(polymers)

monomers = pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/OriginalLibrary_Monomers.csv', header=0)
print(monomers)
with open('Original_Library_Analysis.csv', mode='w') as file:
    file_writer=csv.writer(file, delimiter=',')
    file_writer.writerow(['Polymer_Number','Polymer_Abbr','Biodegradability','Crystallinity'])

for ind in polymers.index:
    number=polymers['Polymer Number'][ind]
    print(number)
    bio=polymers['Biodegradability'][ind]
    print(bio)
    if bio == 'yes':
        biodeg='N'
    if bio=="no":
        biodeg='Y'
    for i in monomers.index:
        if polymers['Monomer 1'][ind]==monomers['Name'][i]:
            n1=monomers['Abbr.'][i]
        elif polymers['Monomer 2'][ind]==monomers['Name'][i]:
            n2=monomers['Abbr.'][i]
    name=str(n1)+'_'+str(n2)
    state=polymers['Physical State'][ind]
    if state=='crystalline':
        cry='Y'
    elif state=='fluid':
        cry='Viscous'
    else:
        cry='N'
    with open('Original_Library_Analysis.csv', mode='a') as file:
        file_writer=csv.writer(file, delimiter=',')
        file_writer.writerow([number,name,biodeg,cry])


