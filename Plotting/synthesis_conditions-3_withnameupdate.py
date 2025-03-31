# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 17:25:07 2023

@author: ChemeGrad2020
"""

###library mass calculations for terpolymers

import matplotlib.pyplot as plt
import numpy as np
import csv
import random
import pandas as pd

#get the data from the random combinations file as a dataframe

combos = pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Library_plan_final.csv')
#drop all the spacing rows
combos=combos.dropna()

#get the data from the monomer info as a dataframe

monomers = pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Monomer Info.csv')
#print(monomers)
#make carboxylic acid dataframe
carbacid=monomers.iloc[:12]
#sort by MW

carbacid_sorted=carbacid.sort_values(by=['MW'])

#make diol dataframe
diol=monomers.iloc[12:]
#sort by bp

diol_sorted=diol.sort_values(by=['Boiling Point (Celcius)- chemspider,760mmHg'])


#start csv file with the synthesis information
with open('Library_Synthesis_final_nameupdate.csv', mode='w') as file:
    file_writer=csv.writer(file, delimiter=',')
    file_writer.writerow(['Polymer_Name_P', 'Monomer', 'Mass_mg', 'mmol', 'equiv', 'MP in C', 'Synthesis Set'])


#decisions for polymerization:
    #1. Set lower mw acid as 500mg (1000mg if only `1)
    #2. Set lower bp diol as 1.1 equiv.
    #3. Select rows to write based on low bp monomers (group by diols) to group by lowest bp

#keep track of which polymers have been set up already    
polymer_numbers=[]

#diol_names = diol_sorted.get(['Carboxylic Acid or Diol'])

#go through all the names in the sorted diols using the index and name column

for ind in diol_sorted.index:
    #check that the monomer is in either of the diol columns
    for i in combos.index:
        ####first check that the combo hasn't already been calculated and move on if it has
        if combos['Polymer Number'][i] in polymer_numbers:
            continue
        if combos['Diol 1'][i]==diol_sorted['Carboxylic Acid or Diol'][ind] or combos['Diol 2'][i] == diol_sorted['Carboxylic Acid or Diol'][ind]:
            #add to the polymer tracking
            polymer_numbers+=[combos['Polymer Number'][i]]
            #### get all the info you need for the stoichiometry questions
            #### the diol with the lowest bp is the monomer you're sorting with
            #generate a dataframe for the carboxylic acids
            #print((combos['Carboxylic Acid 1'][i]))
            CAs=monomers[(monomers['Carboxylic Acid or Diol']==str(combos['Carboxylic Acid 1'][i])) | (monomers['Carboxylic Acid or Diol']==str(combos['Carboxylic Acid 2'][i]))]
            #print(CAs)
            CAs_sorted=CAs.sort_values(by=['MW']) #this makes the first index the lower MW weight one
            #start writing the polymer name
            name=''
            #iterate through all the members of CAs_sorted
            for j in range(len(CAs_sorted.index)):
                if len(CAs_sorted.index)==2:
                    name+=str(CAs_sorted['Polymer Name Abbr.'][CAs_sorted.index[j]])+'25_'
                elif len(CAs_sorted.index)==1:
                    name+=str(CAs_sorted['Polymer Name Abbr.'][CAs_sorted.index[j]])+'50_'
            #set mass for lower MW diacid
            if len(CAs_sorted.index)==2:
                CA_lowmg=500
                #calculate mol basis for reaction
                mmol_base=500/(CAs_sorted['MW'][CAs_sorted.index[0]])
                CA_eql=1
                #calculate mass for high MW diacid
                CA_highmg=mmol_base*CAs_sorted['MW'][CAs_sorted.index[1]]
                CA_eqh=1
            else:
                CA_lowmg=1000
                mmol_base=1000/(CAs_sorted['MW'][CAs_sorted.index[0]])/2
                CA_eql=2
            
            ####switch to diols
            #get dataframe with the relevant diols
            DLs=monomers[(monomers['Carboxylic Acid or Diol']==str(combos['Diol 1'][i])) | (monomers['Carboxylic Acid or Diol']==str(combos['Diol 2'][i]))]
            #sort relative to BP
            DLs_sorted=DLs.sort_values(by=['Boiling Point (Celcius)- chemspider,760mmHg'])
            #the lower bp diol has equivalence of 1.1
            #calculate mass of low bp diol
            if len(CAs_sorted.index)==2:    
                DL_low=mmol_base*DLs_sorted['MW'][DLs_sorted.index[0]]*2.1
                DL_eql=2.1
            else:
                DL_low=mmol_base*DLs_sorted['MW'][DLs_sorted.index[0]]*1.05
                DL_eql=1.05
            #calculate mass of higher bp diol
                DL_high=mmol_base*DLs_sorted['MW'][DLs_sorted.index[1]]
                DL_eqh=1
            #adjust name
            for j in range(len(DLs_sorted.index)):
                if len(DLs_sorted.index)==2:
                    name+=str(DLs_sorted['Polymer Name Abbr.'][DLs_sorted.index[j]])+'25_'
                elif len(DLs_sorted.index)==1:
                    name+=str(DLs_sorted['Polymer Name Abbr.'][DLs_sorted.index[j]])+'50_'
            name+='P'
            #calculate synthesis set
            syn_set=(len(polymer_numbers)//11)+1
            
            #write the information into the csv file
            with open('Library_Synthesis_final_nameupdate.csv', mode='a', newline='') as file:
                file_writer=csv.writer(file, delimiter=',')
                file_writer.writerow([name, CAs_sorted['Carboxylic Acid or Diol'][CAs_sorted.index[0]], CA_lowmg, mmol_base*CA_eql, CA_eql, CAs_sorted['Melting Point (Celcius)-chemspider'][CAs_sorted.index[0]], syn_set])
                if len(CAs_sorted.index)==2:
                    file_writer.writerow([name, CAs_sorted['Carboxylic Acid or Diol'][CAs_sorted.index[1]], CA_highmg, mmol_base*CA_eqh, CA_eqh, CAs_sorted['Melting Point (Celcius)-chemspider'][CAs_sorted.index[1]], syn_set])
                file_writer.writerow([name, DLs_sorted['Carboxylic Acid or Diol'][DLs_sorted.index[0]], DL_low, mmol_base*DL_eql, DL_eql, DLs_sorted['Melting Point (Celcius)-chemspider'][DLs_sorted.index[0]], syn_set])
                if len(DLs_sorted.index)==2:
                    file_writer.writerow([name, DLs_sorted['Carboxylic Acid or Diol'][DLs_sorted.index[1]], DL_high, mmol_base*DL_eqh, DL_eqh, DLs_sorted['Melting Point (Celcius)-chemspider'][DLs_sorted.index[1]], syn_set])
                file_writer.writerow(["","","","","","",""])
                
#check that there are 300 unique polymers
df=pd.read_csv('Library_Synthesis_final_nameupdate.csv',header=0)
df=df.dropna(axis=0,how='all')
uniques=df.drop_duplicates(subset=['Polymer_Name_P'],keep='first')
print(uniques['Polymer_Name_P'])
print(len(uniques.index))