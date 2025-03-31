# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:36:32 2024

@author: ChemeGrad2020
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
import csv
import random
import pandas as pd
import math
mpl.rcParams["font.family"]='Arial'


#repeat unit MW
linears=['SBA','SCA','ADA','AZA','SUA','GLA','OTD','NND','DCD','PPD','PTD','BND','HND','ETG']
linear_numbers=[10,4,6,9,8,5,8,9,10,3,5,4,6,2]

terpolymers = pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Library_Data_Final_crystallinity.csv',header=0)
julia=pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Julia_Data_crystalline.csv', header=0)
monomers = pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Monomer Info_DMTDMI.csv')
original=pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Original_Library_Analysis.csv',header=0)

# print(original)
original_data=[]
terpolymer_data=[]
for i in range(len(original.index)):
    if str(original['Polymer_Abbr'][i])[0:3] in linears and str(original['Polymer_Abbr'][i])[4:7] in linears:
        value=linear_numbers[linears.index(str(original['Polymer_Abbr'][i])[0:3])]+linear_numbers[linears.index(str(original['Polymer_Abbr'][i])[4:7])]
        if original['Biodegradability'].tolist()[i]=='N':
            original_data+=[[value,0]]
        else:
            original_data+=[[value,1]]
#averaging across polymers with only full linear backbones
for i in range(len(terpolymers.index)):
    if str(terpolymers['Polymer_Name_P'][i])[0:3] in linears and str(terpolymers['Polymer_Name_P'][i])[6:9] in linears and str(terpolymers['Polymer_Name_P'][i])[12:15] in linears:
        value=linear_numbers[linears.index(str(terpolymers['Polymer_Name_P'][i])[0:3])]*terpolymers['Ratio 1'][i]+linear_numbers[linears.index(str(terpolymers['Polymer_Name_P'][i])[6:9])]*terpolymers['Ratio 2'][i]+linear_numbers[linears.index(str(terpolymers['Polymer_Name_P'][i])[12:15])]*terpolymers['Ratio 3'][i]
        print(terpolymers["Polymer_Name_P"][i],round(value/2),terpolymers['Biodegradable'][i])
        value=round(value/2)
        if terpolymers['Biodegradable'][i]=='N':
            terpolymer_data+=[[value,0]]
        else:
            terpolymer_data+=[[value,1]]
terpolymer_limited=terpolymer_data
original_array=np.array(original_data)
terpolymer_array_all=np.array(terpolymer_data)
# #taking all polymers with 2 only linear segments
# linear_acids=['SBA','SCA','ADA','AZA','SUA','GLA']
# linear_diols=['OTD','NND','DCD','PPD','PTD','BND','HND','ETG']
# for i in range(len(terpolymers.index)):
#     if str(terpolymers['Polymer_Name_P'][i])[0:3] in linears and str(terpolymers['Polymer_Name_P'][i])[6:9] in linears and str(terpolymers['Polymer_Name_P'][i])[12:15] in linears:
#         value=linear_numbers[linears.index(str(terpolymers['Polymer_Name_P'][i])[0:3])]*terpolymers['Ratio 1'][i]+linear_numbers[linears.index(str(terpolymers['Polymer_Name_P'][i])[6:9])]*terpolymers['Ratio 2'][i]+linear_numbers[linears.index(str(terpolymers['Polymer_Name_P'][i])[12:15])]*terpolymers['Ratio 3'][i]
#         print(terpolymers["Polymer_Name_P"][i],value/2)
#         value=round(value/2)
#         if terpolymers['Biodegradable'][i]=='N':
#             terpolymer_data+=[[value,0]]
#         else:
#             terpolymer_data+=[[value,1]]
#     elif str(terpolymers['Polymer_Name_P'][i])[0:3] in linear_acids and str(terpolymers['Polymer_Name_P'][i])[6:9] in linear_diols:
#         value=linear_numbers[linears.index(str(terpolymers['Polymer_Name_P'][i])[0:3])]+linear_numbers[linears.index(str(terpolymers['Polymer_Name_P'][i])[6:9])]
#         if terpolymers['Biodegradable'][i]=='N':
#             terpolymer_data+=[[value,0]]
#         else:
#             terpolymer_data+=[[value,1]]
#     elif str(terpolymers['Polymer_Name_P'][i])[0:3] in linear_acids and str(terpolymers['Polymer_Name_P'][i])[12:15] in linear_diols:
#         value=linear_numbers[linears.index(str(terpolymers['Polymer_Name_P'][i])[0:3])]+linear_numbers[linears.index(str(terpolymers['Polymer_Name_P'][i])[12:15])]
#         if terpolymers['Biodegradable'][i]=='N':
#             terpolymer_data+=[[value,0]]
#         else:
#             terpolymer_data+=[[value,1]]
#     elif str(terpolymers['Polymer_Name_P'][i])[6:9] in linear_acids and str(terpolymers['Polymer_Name_P'][i])[12:15] in linear_diols:
#         value=linear_numbers[linears.index(str(terpolymers['Polymer_Name_P'][i])[6:9])]+linear_numbers[linears.index(str(terpolymers['Polymer_Name_P'][i])[12:15])]
#         if terpolymers['Biodegradable'][i]=='N':
#             terpolymer_data+=[[value,0]]
#         else:
#             terpolymer_data+=[[value,1]]
original_array=np.array(original_data)
# terpolymer_array_set=np.array(terpolymer_data)
print('original')
print(original_array)
print('terpolymer only full linear')
print(terpolymer_array_all)
# print('terpolymer_all with linear repeat')
# print(terpolymer_array_set)

carbons=np.array([0,8,11,14,17,20])
original_averages=[]

terpolymer_averages_set=[]
terpolymer_avg_total=[]
totals_orig=[]
totals_ter=[]
totals_limit=[]
for m in range(1,len(carbons)):
    j=carbons[m]
    y=carbons[m-1]
    degrades=0
    total=0
    for k in range(len(original_data)):
        if original_data[k][0]<=j and original_data[k][0]>y:
            degrades+=original_data[k][1]
            total+=1
    t_degrades=0
    t_degrades_limit=0
    t_total=0
    t_total_limit=0
    for t in range(len(terpolymer_data)):
        if terpolymer_data[t][0]<=j and terpolymer_data[t][0]>y :
            t_degrades+=terpolymer_data[t][1]
            t_total+=1
    for t in range(len(terpolymer_limited)):
        if terpolymer_limited[t][0]<=j and terpolymer_limited[t][0]>y:
            t_degrades_limit+=terpolymer_limited[t][1]
            t_total_limit+=1
    if total!=0:
        original_averages+=[[j,degrades/total]]
    else:
        original_averages+=[[j,0]]
    # if t_total!=0:
    #     terpolymer_averages_set+=[[j,t_degrades/t_total]]
    # else:
    #     terpolymer_averages_set+=[[j,0]]
    if t_total_limit!=0:
        terpolymer_avg_total+=[[j,t_degrades_limit/t_total_limit]]
    else:
        terpolymer_avg_total+=[[j,0]]
    totals_orig+=[total]
    # totals_ter+=[t_total]
    totals_limit+=[t_total_limit]
original_averages_np=np.array(original_averages)
# terpolymer_averages_np=np.array(terpolymer_averages_set)
terpolymer_total_np=np.array(terpolymer_avg_total)
# print(terpolymer_total_np)
print(terpolymer_limited)
#print(totals_orig)

X = ['< 8','9-11','12-14','15-17','18-20'] 
# Y_ter = terpolymer_averages_np[:,1].tolist()
Y_tot =terpolymer_total_np[:,1].tolist()
Z_original = original_averages_np[:,1].tolist()
  
X_axis = np.arange(len(X)) 


plt.figure()      
plt.bar(X_axis - 0.2, Y_tot, 0.4, label = 'Copolymers',color='#527b08')
for i in X_axis:
        plt.text(i-0.2, Y_tot[i], totals_limit[i], ha = 'center',color='#527b08',size=14)
plt.bar(X_axis + 0.2, Z_original, 0.4, label = 'Reduced Original Library',color='#debd29')
for i in X_axis:
        plt.text(i+0.2, Z_original[i], totals_orig[i], ha= 'center',color='black',size=14) 
plt.xticks(X_axis, X,size=14) 
plt.yticks(size=14)
plt.xlabel("Carbon Length of Linear Repeat Unit", size=14) 
plt.ylabel("Fraction of Polymers\nwhich Biodegrade",size=14) 

plt.legend(fontsize=14) 

plt.savefig('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/linear_repeat_unit_only_linear_fixed_colors.pdf',dpi=600, bbox_inches='tight')     
plt.show()



######################## toxicity stuff
bpa_total=0
bpa_no_growth=0
for i in range(len(terpolymers.index)):
    if 'BPA' in str(terpolymers['Polymer_Name_P'][i]):
        bpa_total+=1
        if terpolymers['Degradation Type'][i]=='No Growth':
            bpa_no_growth+=1
for i in range(len(julia.index)):
    if 'BPA' in str(julia['Polymer_Name_P'][i]):
        bpa_total+=1
        if julia['Degradation Type'][i]=='No Growth':
            bpa_no_growth+=1
print(bpa_no_growth/bpa_total)

####################### looking at polymer degradability compared to terpolymer
deg_ter_bothY=0
deg_ter_oneYoneN=0
deg_ter_twoN=0
nondeg_ter_bothY=0
nondeg_ter_oneYoneN=0
nondeg_ter_twoN=0
no_data_biodeg=0
no_data_nonbiodeg=0

nondeg_nogrowth_Y=0
nondeg_growth_Y=0
nondeg_nogrowth_N=0
nondeg_growth_N=0
nondeg_nogrowth_YY=0
nondeg_nogrowth_YN=0
nondeg_nogrowth_NN=0
nondeg_nogrowth_ND=0

for i in range(len(terpolymers.index)):
    mon1=str(terpolymers['Polymer_Name_P'][i])[0:3]
    mon2=str(terpolymers['Polymer_Name_P'][i])[6:9]
    mon3=str(terpolymers['Polymer_Name_P'][i])[12:15]
    value1=None
    value2=None
    for j in range(len(original.index)):
        if mon1 in str(original['Polymer_Abbr'][j]) and mon3 in str(original['Polymer_Abbr'][j]):
            value1=original['Biodegradability'][j]
        if mon1 in str(original['Polymer_Abbr'][j]) and mon2 in str(original['Polymer_Abbr'][j]):
            value2=original['Biodegradability'][j]
        elif mon2 in str(original['Polymer_Abbr'][j]) and mon3 in str(original['Polymer_Abbr'][j]):
            value2=original['Biodegradability'][j]
    if value1==None or value2==None:
        if terpolymers['Biodegradable'][i]=='Y':
            no_data_biodeg+=1
        else:
            no_data_nonbiodeg+=1
            if terpolymers['Degradation Type'][i]=='No Growth' or terpolymers['Degradation Type'][i]=='Hydrolysis':
                nondeg_nogrowth_ND+=1
    elif value1==None and value2==None:
        if terpolymers['Biodegradable'][i]=='Y':
            no_data_biodeg+=1
        else:
            no_data_nonbiodeg+=1
    elif terpolymers['Biodegradable'][i]=='Y':
        if value1=='Y' and value2=='Y':
            deg_ter_bothY+=1
        elif (value1=='Y' and value2=='N'):
            deg_ter_oneYoneN+=1
        elif (value1=='N' and value2=='Y'):
            deg_ter_oneYoneN+=1
        elif value1=='N' and value2=='N':
            deg_ter_twoN+=1
    elif terpolymers['Biodegradable'][i]=='N':
        if value1=='Y' and value2=='Y':
            nondeg_ter_bothY+=1
            if terpolymers['Degradation Type'][i]=='No Growth' or terpolymers['Degradation Type'][i]=='Hydrolysis':
                nondeg_nogrowth_Y+=1
                nondeg_nogrowth_YY+=1
            else:
                nondeg_growth_Y+=1
            print(terpolymers['Polymer_Name_P'][i],terpolymers['Degradation Type'][i],value1,value2)
        elif (value1=='Y' and value2=='N'):
            nondeg_ter_oneYoneN+=1
            print('one',terpolymers['Polymer_Name_P'][i],terpolymers['Degradation Type'][i],value1,value2)
            if terpolymers['Degradation Type'][i]=='No Growth' or terpolymers['Degradation Type'][i]=='Hydrolysis':
                nondeg_nogrowth_Y+=1
                nondeg_nogrowth_YN+=1
            else:
                nondeg_growth_Y+=1
        elif (value1=='N' and value2=='Y'):
            nondeg_ter_oneYoneN+=1
            print('one',terpolymers['Polymer_Name_P'][i],terpolymers['Degradation Type'][i],value1,value2)
            if terpolymers['Degradation Type'][i]=='No Growth' or terpolymers['Degradation Type'][i]=='Hydrolysis':
                nondeg_nogrowth_Y+=1
                nondeg_nogrowth_YN+=1
            else:
                nondeg_growth_Y+=1
        elif value1=='N' and value2=='N':
            nondeg_ter_twoN+=1
            if terpolymers['Degradation Type'][i]=='No Growth' or terpolymers['Degradation Type'][i]=='Hydrolysis':
                nondeg_nogrowth_N+=1
                nondeg_nogrowth_NN+=1
            else:
                nondeg_growth_N+=1
for i in range(len(julia.index)):
    mon1=str(julia['Polymer_Name_P'][i])[0:3]
    mon2=str(julia['Polymer_Name_P'][i])[6:9]
    mon3=str(julia['Polymer_Name_P'][i])[12:15]
    value1=None
    value2=None
    for j in range(len(original.index)):
        if mon1 in str(original['Polymer_Abbr'][j]) and mon3 in str(original['Polymer_Abbr'][j]):
            value1=original['Biodegradability'][j]
        if mon1 in str(original['Polymer_Abbr'][j]) and mon2 in str(original['Polymer_Abbr'][j]):
            value2=original['Biodegradability'][j]
        elif mon2 in str(original['Polymer_Abbr'][j]) and mon3 in str(original['Polymer_Abbr'][j]):
            value2=original['Biodegradability'][j]
    if value1==None or value2==None:
        if terpolymers['Biodegradable'][i]=='Y':
            no_data_biodeg+=1
        else:
            no_data_nonbiodeg+=1
    elif value1==None and value2==None:
        if terpolymers['Biodegradable'][i]=='Y':
            no_data_biodeg+=1
        else:
            no_data_nonbiodeg+=1
    elif julia['Biodegradable'][i]=='Y':
        if value1=='Y' and value2=='Y':
            deg_ter_bothY+=1
        elif (value1=='Y' and value2=='N'):
            deg_ter_oneYoneN+=1
        elif (value1=='N' and value2=='Y'):
            deg_ter_oneYoneN+=1
        elif value1=='N' and value2=='N':
            deg_ter_twoN+=1
    elif julia['Biodegradable'][i]=='N':
        if value1=='Y' and value2=='Y':
            nondeg_ter_bothY+=1
            print(julia['Polymer_Name_P'][i])
            if julia['Degradation Type'][i]=='No Growth' or julia['Degradation Type'][i]=='Hydrolysis':
                nondeg_nogrowth_Y+=1
                nondeg_nogrowth_YY+=1
            else:
                nondeg_growth_Y+=1
        elif (value1=='Y' and value2=='N'):
            nondeg_ter_oneYoneN+=1
            print('one',julia['Polymer_Name_P'][i],julia['Degradation Type'][i],value1,value2)
            if julia['Degradation Type'][i]=='No Growth' or julia['Degradation Type'][i]=='Hydrolysis':
                nondeg_nogrowth_Y+=1
                nondeg_nogrowth_YN+=1
            else:
                nondeg_growth_Y+=1
        elif (value1=='N' and value2=='Y'):
            nondeg_ter_oneYoneN+=1
            print('one',julia['Polymer_Name_P'][i],julia['Degradation Type'][i],value1,value2)
            if julia['Degradation Type'][i]=='No Growth' or julia['Degradation Type'][i]=='Hydrolysis':
                nondeg_nogrowth_Y+=1
                nondeg_nogrowth_YN+=1
            else:
                nondeg_growth_Y+=1
        elif value1=='N' and value2=='N':
            nondeg_ter_twoN+=1
            if julia['Degradation Type'][i]=='No Growth' or julia['Degradation Type'][i]=='Hydrolysis':
                nondeg_nogrowth_N+=1
                nondeg_nogrowth_NN+=1
            else:
                nondeg_growth_N+=1
two_bi_deg=[deg_ter_bothY,nondeg_ter_bothY]
one_one=[deg_ter_oneYoneN,nondeg_ter_oneYoneN]
two_N=[deg_ter_twoN,nondeg_ter_twoN]
no_data=[no_data_biodeg,no_data_nonbiodeg]

print(two_bi_deg,one_one,two_N,no_data_biodeg, no_data_nonbiodeg)

print('deg_ter_both Y', (deg_ter_bothY))
print('deg_ter_YN',deg_ter_oneYoneN)
print('deg_ter_NN',deg_ter_twoN)
print('nondeg_ter_YY',(nondeg_ter_bothY))
print('nondeg_ter_YN',(nondeg_ter_oneYoneN))
print('nondeg_ter_NN', nondeg_ter_twoN)
print('no_data (deg,nondeg)',no_data)
colors=['#debd29','#8bbd31','#527b08','#295208']
fig,ax=plt.subplots()
ax.barh(np.array([1,2]),np.array(two_bi_deg),facecolor=colors[3],edgecolor=colors[3],label='Two')
# ax.barh([2],nondeg_nogrowth_YY,facecolor=colors[3],edgecolor='darkgrey',hatch=r"//",label='No Colony Growth')
ax.barh(np.array([1,2]),np.array(one_one),left=np.array(two_bi_deg),facecolor=colors[1],edgecolor=colors[1],label='One')
ax.barh([2],nondeg_nogrowth_YN,left=two_bi_deg[1],facecolor=colors[1],edgecolor='darkgrey',hatch=r"//")
ax.barh(np.array([1,2]),np.array(two_N),left=np.array(two_bi_deg)+np.array(one_one),facecolor=colors[0],edgecolor=colors[0],label='Zero')
ax.barh([2],nondeg_nogrowth_NN,left=two_bi_deg[1]+one_one[1],facecolor=colors[0],edgecolor='darkgrey',hatch=r"//")
ax.barh(np.array([1,2]),np.array(no_data),left=np.array(two_bi_deg)+np.array(one_one)+np.array(two_N),facecolor='grey',edgecolor='grey',label='No Data')
ax.barh([2],nondeg_nogrowth_YY,facecolor='white',edgecolor='darkgrey',hatch=r"//",label='No Colony Growth')
ax.barh(([2]),nondeg_nogrowth_ND,left=two_bi_deg[1]+one_one[1]+two_N[1],facecolor='grey',edgecolor='darkgrey',hatch=r"//")
ax.set_yticks([1,2],['Biodegradable\nTerpolymers','Non-Biodegradable\nTerpolymers'],fontsize=14)
ax.set_xticks([0,50,100,150,200,250,300],fontsize=14)
plt.xticks(size=14)
ax.set_xlabel('Number of Terpolymers',fontsize=14)
ax.legend(title='Number of Terpolymer Repeat Units\nwhich Biodegrade as Binary Copolymers',ncols=2, fontsize=14,title_fontsize=14)
fig.show()
# plt.savefig('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/biodeg_repeat_unit_breakdown_withnogrowth.pdf',dpi=600, bbox_inches='tight')
print('nondegnogrowthyn',nondeg_nogrowth_YN)
#look at degradation behavior of the non-degrading polymers
nondeg_Ytotal=nondeg_nogrowth_Y+nondeg_growth_Y
nondeg_Ntotal=nondeg_nogrowth_N+nondeg_growth_N
fig,ax=plt.subplots()
print('total nondeg',(nondeg_Ytotal+nondeg_Ntotal))
print('colony behavior nondeg polymers',nondeg_nogrowth_Y/nondeg_Ytotal, nondeg_nogrowth_N/nondeg_Ntotal)
ax.bar([1,2],[nondeg_nogrowth_Y/nondeg_Ytotal,nondeg_nogrowth_N/nondeg_Ntotal],color=colors[0],label='No colony growth')
ax.bar([1,2],[nondeg_growth_Y/nondeg_Ytotal,nondeg_growth_N/nondeg_Ntotal],bottom=[nondeg_nogrowth_Y/nondeg_Ytotal,nondeg_nogrowth_N/nondeg_Ntotal],color=colors[3],label='Colony growth')  
ax.set_xticks([1,2],['Biodegradable consituent\nbinary copolymer(s)','No biodegradable\nconstituent binary copolymers'])                                                                                           
plt.xticks(size=14)
ax.set_ylabel('Fraction of Terpolymers',rotation=90,fontsize=14)
ax.set_yticks([0,0.2,0.4,0.6,0.8,1],fontsize=14)
ax.legend(fontsize=14)
fig.show()
ax.tick_params(axis='y',labelsize=14)
# plt.savefig('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/nogrowth_non_degrading_terpolymers.pdf',dpi=600, bbox_inches='tight')
num_PTD_polymers=0
num_PTD_nogrowth=0
num_bpa_polymers=0
num_bpa_nogrowth=0
total_nogrowth=0
num_gge_polymers=0
num_gge_nogrowth=0
num_DMT_polymers=0
num_DMT_BPA_nogrowth=0
num_DMT_GGE_nogrowth=0
num_DMT_other_nogrowth=0
num_DMI_polymers=0
num_DMI_BPA_nogrowth=0
num_DMI_GGE_nogrowth=0
num_DMI_other_nogrowth=0
for i in range(len(terpolymers.index)):
    if terpolymers['Degradation Type'][i]=='No Growth' or terpolymers['Degradation Type'][i]=='Hydrolysis':
        total_nogrowth+=1
        print('no growth',terpolymers['Polymer_Name_P'][i])
    if 'BPA' in str(terpolymers['Polymer_Name_P'][i]):
        num_bpa_polymers+=1
        if terpolymers['Degradation Type'][i]=='No Growth' or terpolymers['Degradation Type'][i]=='Hydrolysis':
            num_bpa_nogrowth+=1
    if 'GGE' in str(terpolymers['Polymer_Name_P'][i]):
        num_gge_polymers+=1
        if terpolymers['Degradation Type'][i]=='No Growth' or terpolymers['Degradation Type'][i]=='Hydrolysis':
            num_gge_nogrowth+=1
    if 'PTD' in str(terpolymers['Polymer_Name_P'][i]):
        num_PTD_polymers+=1
        if terpolymers['Degradation Type'][i]=='No Growth' or terpolymers['Degradation Type'][i]=='Hydrolysis':
            num_PTD_nogrowth+=1
print('bpa',num_bpa_polymers)
print('gge',num_gge_polymers)
print(len(julia.index))      
for i in range(len(julia.index)):
    if julia['Degradation Type'][i]=='No Growth' or julia['Degradation Type'][i]=='Hydrolysis':
        total_nogrowth+=1
        print('no growth',julia['Polymer_Name_P'][i])
    if 'BPA' in str(julia['Polymer_Name_P'][i]):
        num_bpa_polymers+=1
    elif 'GGE' in str(julia['Polymer_Name_P'][i]):
        num_gge_polymers+=1
    elif 'PTD' in str(julia['Polymer_Name_P'][i]):
        num_PTD_polymers+=1
    if 'DMT' in str(julia['Polymer_Name_P'][i]):
        num_DMT_polymers+=1
        if julia['Degradation Type'][i]=='No Growth' or julia['Degradation Type'][i]=='Hydrolysis':
            if 'BPA' in julia['Polymer_Name_P'][i]:
                num_DMT_BPA_nogrowth+=1
            elif 'GGE' in julia['Polymer_Name_P'][i]:
                num_DMT_GGE_nogrowth+=1
            elif 'PTD' in julia['Polymer_Name_P'][i]:
                num_PTD_nogrowth+=1
            else:
                num_DMT_other_nogrowth+=1
    if 'DMI' in str(julia['Polymer_Name_P'][i]):
        num_DMI_polymers+=1
        if julia['Degradation Type'][i]=='No Growth' or julia['Degradation Type'][i]=='Hydrolysis':
            if 'BPA' in julia['Polymer_Name_P'][i]:
                num_DMI_BPA_nogrowth+=1
            elif 'GGE' in julia['Polymer_Name_P'][i]:
                num_DMI_GGE_nogrowth+=1
            else:
                num_DMI_other_nogrowth+=1
print('total',total_nogrowth)
fig,ax=plt.subplots()
plt.rcParams["hatch.linewidth"] = 4
ax.tick_params(axis='x',labelsize=14)
print((num_bpa_nogrowth+num_DMT_BPA_nogrowth+num_DMI_BPA_nogrowth+num_DMI_other_nogrowth+num_gge_nogrowth))
ax.barh(1,(num_bpa_nogrowth)/total_nogrowth,height=0.4,facecolor=colors[3],edgecolor=colors[3],label="BPA")
ax.barh(1,num_DMT_BPA_nogrowth/total_nogrowth,height=0.4,left=(num_bpa_nogrowth)/total_nogrowth,facecolor=colors[3],edgecolor=colors[2],hatch=r"//",label="BPA and Terephthalate")
ax.barh(1,num_DMI_BPA_nogrowth/total_nogrowth,height=0.4,left=(num_bpa_nogrowth+num_DMT_BPA_nogrowth)/total_nogrowth,facecolor=colors[3],edgecolor=colors[1],hatch=r"//",label="BPA and Isophthalate")
ax.barh(1,num_DMI_other_nogrowth/total_nogrowth,height=0.4,left=(num_bpa_nogrowth+num_DMT_BPA_nogrowth+num_DMI_BPA_nogrowth)/total_nogrowth,facecolor=colors[1],edgecolor=colors[1],label="Isophthalate")
ax.barh(1,num_gge_nogrowth/total_nogrowth,height=0.4,left=(num_bpa_nogrowth+num_DMT_BPA_nogrowth+num_DMI_BPA_nogrowth+num_DMI_other_nogrowth)/total_nogrowth,facecolor=colors[0],edgecolor=colors[0],label="Guaicol Glycerol Ether")
ax.barh(1,num_PTD_nogrowth/total_nogrowth,height=0.4,left=(num_bpa_nogrowth+num_DMT_BPA_nogrowth+num_DMI_BPA_nogrowth+num_DMI_other_nogrowth+num_gge_nogrowth)/total_nogrowth,facecolor='grey',edgecolor='grey',label="1,5-Pentanediol")
print(total_nogrowth)
ax.legend(fontsize=14,loc='upper center',bbox_to_anchor=(0.5,-0.15),ncol=2)
ax.set_xticks([0,0.2,0.4,0.6,0.8,1],fontsize=14)
ax.set_xlim([0,1])
ax.set_xlabel("Fraction of Polymers with No Colony Growth",fontsize=14)
ax.set_yticks([])
fig.show()
plt.savefig('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/nogrowth_aromatic_terpolymers_legend_sep.pdf',dpi=600, bbox_inches='tight')




#####looking at crystallinity


crystalline=0
crystalline_deg=0
amorphous=0
amorphous_deg=0
solid=0
solid_deg=0
for j in range(len(terpolymers.index)):
    i=str(terpolymers['Crystalline'][j])
    print(i)
    if i == 'Y':
        crystalline+=1
        if terpolymers['Biodegradable'][j]=='Y':
            crystalline_deg+=1
    if i =='N':
        solid+=1
        if terpolymers['Biodegradable'][j]=='Y':
            solid_deg+=1
    if i =='nan':
        amorphous+=1
        if terpolymers['Biodegradable'][j]=='Y':
            amorphous_deg+=1
for j in range(len(julia.index)):
    i=str(julia['Crystalline'][j])
    if i == 'Y':
        crystalline+=1
        if julia['Biodegradable'][j]=='Y':
            crystalline_deg+=1
    if i =='N':
        solid+=1
        if julia['Biodegradable'][j]=='Y':
            solid_deg+=1
    if i =='nan':
        amorphous+=1
        if julia['Biodegradable'][j]=='Y':
            amorphous_deg+=1
ogcrystalline=0
ogamorphous=0
ogsolid=0      
for i in original['Crystallinity']:
    # print(i)
    if i == 'Y':
        ogcrystalline+=1
    if i =='N':
        ogsolid+=1
    if i =='Viscous':
        ogamorphous+=1
print(crystalline_deg/crystalline,amorphous_deg/amorphous,solid_deg/solid,crystalline)
# print(ogcrystalline/(ogamorphous+ogsolid+ogcrystalline))

fig,ax=plt.subplots()
ax.barh([0,0.5,1],[crystalline,solid,amorphous],height=0.4,facecolor=colors[3],edgecolor=colors[3],label='Biodegradable')
ax.barh([0,0.5,1],[crystalline-crystalline_deg,solid-solid_deg,amorphous-amorphous_deg],left=[crystalline_deg,solid_deg,amorphous_deg],height=0.4,facecolor=colors[3],edgecolor=colors[1],hatch=r"//",label='Non-Biodegradable')
ax.legend(fontsize=14,loc='lower right')
print(crystalline_deg/crystalline,solid_deg/solid,amorphous_deg/amorphous)
print('polymer state',crystalline, solid, amorphous)
ax.set_xticks([0,25,50,75,100,125,150,175],fontsize=14)
ax.set_yticks([0,0.5,1],['Crystalline','Non-Crystalline\nSolid','Viscous'],fontsize=14)
ax.tick_params(axis='both',labelsize=14)
ax.set_xlabel('Number of Polymers',fontsize=14)
fig.show()
plt.savefig('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/biodegradability_crystallinity_polymernum.png',dpi=600, bbox_inches='tight')

####### count number of polymers with sulfur
####### count number of polymers with oxygen
oxygen=['DEG','TEG','EBA','DGA','TPG','DPG']
sulfur=['TOA']
oxy=0
oxyar=0
for mon in sulfur:
    for name in terpolymers['Polymer_Name_P']:
        if mon in name:
            oxy+=1
    for name in julia['Polymer_Name_P']:
        if mon in name:
            oxyar+=1
print(oxy,oxyar)


#terephthalic isophthalic odd even

list_interest=['ETG','PPD','BND','PTD','HND','HEP','OTD','NND','DCD']
#make carboxylic acid dataframe
carbacid=monomers.iloc[:12]
print(carbacid)
# carbacid = carbacid.drop([0,2])
carbacid=carbacid.sort_values(by=['MW'])
print(carbacid)
#make diol dataframe
diol=monomers.iloc[12:36]
diol=diol.sort_values(by=['MW'])
print(diol)


listnp=[]
labels=[]
acids_list = carbacid['Polymer Name Abbr.'].tolist()+['HEP']
# acids_list=acids_list1[0:9]+['SBA']
print(acids_list)
diols_list=diol['Polymer Name Abbr.'].tolist()
print(diols_list)
combined=acids_list+diols_list
# combined=['DMT','DMI']+diols_list
acids_names=carbacid['SMILES'].tolist()+['OCCCCCCCO']
# acids_names=acids_names[0:9]+['Sebacic Acid']
diols_names=diol['SMILES'].tolist()
combined_names=acids_names+diols_names
print(combined_names)
smiles_interest=[]
for mon in list_interest:
    smiles_interest+=[combined_names[combined.index(mon)]]
    print(mon,[combined_names[combined.index(mon)]] )

deg=np.zeros(len(list_interest))
tot=np.zeros(len(list_interest))
cry=np.zeros(len(list_interest))
both=np.zeros(len(list_interest))
for mon in range(len(list_interest)):
    for j in range(len(julia.index)):
        if list_interest[mon] in julia['Polymer_Name_P'][j]:
            tot[mon]+=1
            if julia['Biodegradable'][j]=='Y':
                deg[mon]+=1
            if julia['Crystalline'][j]=='Y':
                cry[mon]+=1
            if julia['Biodegradable'][j]=='Y' and julia['Crystalline'][j]=='Y':
                both[mon]+=1
print(tot,deg,cry,both)

fig, ax= plt.subplots()
xcord=np.arange(0,len(list_interest))
ax.bar(xcord,deg/tot,color=colors[3])
ax.set_xticks(xcord,smiles_interest,fontsize=14)
plt.setp(ax.get_xticklabels(), rotation=30, ha="right",
          rotation_mode="anchor")
ax.set_yticks([0.2,0.4,0.6,0.8,1],fontsize=14)
ax.tick_params(axis='y',labelsize=14)
ax.set_ylabel('Fraction of Polymers\nwhich Biodegrade',fontsize=14)
plt.savefig('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/phthalate_linear_degradation.pdf',dpi=600, bbox_inches='tight')
fig.show()
