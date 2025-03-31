# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 16:23:10 2023

@author: kfransen
"""

###library design for terpolymers

import matplotlib.pyplot as plt
import numpy as np
import csv
import random
import pandas as pd

###potential carboxylic acid monomer options
carbacid=['terephthalic acid', 'diglycolic acid', 'sebacic acid', 'succinic acid'
          , 'isophthalic acid', 'adipic acid', 'azelaic acid','suberic acid', '2,2-Thiodiacetic Acid'
          , '1,4-Cyclohexanedicarboxylic Acid','Glutaric Acid','2,2-[ethylenebis(oxy)] bisacetic acid']

diols= ['bisphenol A', '2,2,4-Trimethyl-1,3-pentanediol'
        , '2,2-Dimethyl-1,3-propanediol', '2,3-Butanediol', '1,8-octanediol'
        , '1,9-nonanediol','1,10-decanediol', '1,4-cyclohexanediol'
        , 'diethylene glycol', 'triethylene glycol', 'Guaiacol glyceryl ether'
        , '2,2-Diethyl-1,3-propanediol', '3-methyl-1,5-pentanediol',
        'tripropylene glycol','1,3-propanediol', '1,5-pentanediol'
        , 'dipropylene glycol','1,4-butanediol','1,6-hexanediol'
        , '2-Methyl-1,3-propanediol','2-Butyl-2-ethyl-1,3-propanediol'
        , 'Ethylene Glycol', '1,2-propanediol'
        ,'1,4-Cyclohexanedimethanol']
#print(len(carbacid)) #there are 12

#print(len(diols)) #there are 24
names=carbacid +diols +['None']
names_nonone=carbacid+diols
print((names))
print(len(names))
#read in monomer file
monomers = pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Monomer Info.csv')
#make carboxylic acid dataframe
carbacids=monomers.iloc[:12]
#sort by MW

carbacid_sorted=carbacids.sort_values(by=['MW'])

#make diol dataframe
diol=monomers.iloc[12:]
#sort by bp

diol_sorted=diol.sort_values(by=['Boiling Point (Celcius)- chemspider,760mmHg'])



#choose a number of polymers here:

polyms=300

# #open a  figure for plotting
# fig, ax=plt.subplots()
# #for the xaxis: random number from 0 to 11 for the carboxylic acids
# #for the xaxis: random number from 12 to 35
# #generate the range of x values and associate them with the labels
# #make range of values
# x=np.arange(0,36,1)
# #print(x)
# #associate with labels
# plt.xticks(x,names_nonone)
# #move the labels to the top
# ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
# #tilt the labels so they're readable
# plt.xticks(rotation=90)
# fig.set_figheight(20)
# fig.set_figwidth(6)
# #write the column names and start file
# with open('Library_plan.csv', mode='w') as file:
#     file_writer=csv.writer(file, delimiter=',')
#     file_writer.writerow(['Polymer Number','Carboxylic Acid 1','Carboxylic Acid 2','Diol 1','Diol 2','Name'])
# #generate all the polymer sets, plot them, and add to csv file
# #start your polymer combo counter at 0
# n=0
# headers={'Name':[]}
# track=pd.DataFrame(headers)
# #print(track)
# #print(track)
# while n< polyms:
#     #generate random numbers for two carboxylic acids
#     one=random.randint(0,11)
#     #print(one)
#     #this has you pick two carboxylic acids for polymers 201-300
#     if n>=200:
#         two=random.randint(0,11)
#         #ensure that carboxylic acids aren't the same
#         while two==one:
#             two=random.randint(0,11)
#     #generate random numbers for two diols
#     else:
#         two=36
#     three=random.randint(12,35)
#     #this has you pick two diols for polymers 1-200
#     if n<200:
#         four=random.randint(12,35)
#         #ensure that diols aren't the same
#         while four==three:
#             four=random.randint(12,35)
#     else:
#         four=36
#     #check that the polymer is unique
#     #first write the name abbr. of the polymer
#     #get the carboxylic acid dataframe
#     CAs=monomers[(monomers['Carboxylic Acid or Diol']==names[one]) | (monomers['Carboxylic Acid or Diol']==names[two])]
#     #sort relative to MW
#     CAs_sorted=CAs.sort_values(by=['MW'])
#     #start writing the name for this file
#     namep=''
#     #iterate through all the members of CAs_sorted
#     for j in range(len(CAs_sorted.index)):
#         namep+=str(CAs_sorted['Polymer Name Abbr.'][CAs_sorted.index[j]])+'_'
#     #get the diol info, sort, and add to the polymer name
#     DLs=monomers[(monomers['Carboxylic Acid or Diol']==names[three]) | (monomers['Carboxylic Acid or Diol']==names[four])]
#     #sort relative to BP
#     DLs_sorted=DLs.sort_values(by=['Boiling Point (Celcius)- chemspider,760mmHg'])
#     for j in range(len(DLs_sorted.index)):
#         namep+=str(DLs_sorted['Polymer Name Abbr.'][DLs_sorted.index[j]])+'_'
#     namep+='P'
#     #add the new row to the tracking dataframe
#     track=track.append({'Name':namep},ignore_index=True)
#     #check for duplicates and drop duplicates keeping only the first
#     track=track.drop_duplicates(subset='Name',keep='first')
#     #calculate the index of the dataframe
#     num=len(track.index)
#     if num==n+1:
#         n=num
#         if n<=200:
#             plt.scatter([one],[n],c='blue',marker='o')
#             plt.scatter([three,four],[n,n],c='red',marker='+')
#         if n>200:
#             plt.scatter([one,two],[n,n],c='blue',marker='o')
#             plt.scatter([three],[n],c='red',marker='+')
#         with open('library_plan.csv', mode='a') as file:
#             file_writer=csv.writer(file,delimiter=',')
#             file_writer.writerow([str(n),names[one],names[two],names[three],names[four],namep])
         
#analysis of the generated library
#import into pandas
df=pd.read_csv('Library_plan_final.csv', header=0)


carboxylicacids1=df.groupby(['Carboxylic Acid 1']).size()
carboxylicacids2=df.groupby([ 'Carboxylic Acid 2']).size()
carboxylicacids=carboxylicacids1+carboxylicacids2
diolsoccur1=df.groupby(['Diol 1']).size()
diolsoccur2=df.groupby(['Diol 2']).size()
dioloccur=diolsoccur1+diolsoccur2

print(carboxylicacids)
print(dioloccur) 






