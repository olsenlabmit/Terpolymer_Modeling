# -*- coding: utf-8 -*-
"""
Created on Wed May 29 18:54:40 2024

@author: kfransen
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
import csv
import random
import pandas as pd
import math
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import StrMethodFormatter
from heatmap import heatmap
from heatmap import annotate_heatmap

mpl.rcParams["font.family"]='Arial'


terpolymers = pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Library_Data_Final_Ratios.csv', header=0)
julia=pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Library_Data_Julia_Ratios.csv', header=0)
monomers = pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Monomer Info_DMTDMI.csv')
original=pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Original_Library_Analysis.csv',header=0)

RF_original_incr=pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/RF_original_training_incr_vs_retrain.csv',header=0)
RF_reduced_incr=pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/RF_reduced_training_incr_vs_retrain.csv',header=0)
SDG_original_incr=pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/SDG_original_training_incr_vs_retrain.csv',header=0)
SDG_reduced_incr=pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/SDG_reduced_training_incr_vs_retrain.csv',header=0)

data_frame_list=[RF_original_incr,RF_reduced_incr,SDG_original_incr,SDG_reduced_incr]
additional_datapoint=[0,25,23,25,24,25,25,28,27,30,31]
total_datapoints=np.array([0,25,48,73,97,122,147,175,202,232,263])
fraction_original=total_datapoints/(total_datapoints + 546*np.ones(11))
fraction_reduced=total_datapoints/(total_datapoints + 197*np.ones(11))

colors=['#debd29','#8bbd31','#527b08','#295208']
markers=['v','^','d']

original_SDG_cond=[['hinge',0.1],['log_loss',0.01],['perceptron',0.001]]
reduced_SDG_cond=[['hinge',0.1],['log_loss',0.1],['perceptron',0.0001]]
original_RF_cond=[[8,4,8],[16,16,128],[16,128,512],[32,256,8]]
reduced_RF_cond=[[8,2,8],[16,4,128],[64,256,8],[256,4,256]]

all_conds=[original_RF_cond,reduced_RF_cond,original_SDG_cond,reduced_SDG_cond]
all_conds_names=['RF_orig','RF_Reduc','SDG_orig','SDG_red']
xvals=[fraction_original,fraction_reduced,fraction_original,fraction_reduced]
#plotting the retrained models
for k in range(len(all_conds)):
    for j in range(len(all_conds[k])):
        dataframe=data_frame_list[k]
        conditions=all_conds[k][j]
        new_df=dataframe.loc[dataframe['Original Model Parameters']==str(conditions)]
        new_df=new_df.loc[new_df['Retrain or Addition'].isin(["Retrain",'Basis'])]
        if k <2:
            new_df=new_df.loc[new_df['Total Trees']==conditions[0]]
            title='Random forest model'
        else:
            title='Linear model'
        fig, ax=plt.subplots(layout='constrained')
        ax.scatter(xvals[k],new_df['Binary Polymers Score'],c=colors[0],marker=markers[0],label='Homopolymers')
        ax.scatter(xvals[k],new_df['Terpolymer Test Score'],c=colors[1],marker=markers[1],label='Copolymers')
        ax.scatter(xvals[k],new_df['Combined Test Score'],c=colors[2],marker=markers[2],label='Combined Polymers')
        plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places
        plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places
        ax.set_xlabel('Fraction of Copolymers in Training Data',fontsize=11)
        ax.set_ylabel('Test Set Accuracy',fontsize=11)
        ax.set_ylim([0.5,1.0])
        ax.set_yticks([0.5,0.6,0.7,0.8,0.9,1.0])
        ax.tick_params(axis='both',which='major',labelsize=11)
        ax.legend(loc='best',ncols=1,fontsize=11)
        ax.tick_params(axis='both',labelsize=11)
        ax.set_title(title,fontsize=11)
        # fig.set_size_inches(3, 3)
        plt.savefig('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/Retrain_'+all_conds_names[k]+'_'+str(conditions)+'.pdf',dpi=300, bbox_inches='tight')
#plotting the models with the additions        
for k in range(len(all_conds)):
    for j in range(len(all_conds[k])):
        dataframe=data_frame_list[k]
        conditions=all_conds[k][j]
        print(conditions)
        new_df=dataframe.loc[dataframe['Original Model Parameters']==str(conditions)]
        new_df=new_df.loc[new_df['Retrain or Addition'].isin(["Addition",'Basis'])]
        fig, ax=plt.subplots(layout='constrained')
        ax.scatter(xvals[k],new_df['Binary Polymers Score'],c=colors[0],marker=markers[0],label='Binary Copolymers')
        ax.scatter(xvals[k],new_df['Terpolymer Test Score'],c=colors[1],marker=markers[1],label='Terpolymers')
        ax.scatter(xvals[k],new_df['Combined Test Score'],c=colors[2],marker=markers[2],label='Combined Polymers')
        plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places
        plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places
        ax.set_xlabel('Fraction of Terpolymers in Training Data')
        ax.set_ylabel('Test Set Accuracy')
        ax.legend(loc='best',ncols=1,fontsize='small')
        # plt.savefig('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/Addition_'+all_conds_names[k]+'_'+str(conditions)+'.png',dpi=300, bbox_inches='tight')       

#plotting accuracy vs landscape
landscapes=pd.read_csv('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/Landscape_study.csv',header=0)
original=landscapes.loc[landscapes['Landscape Category']=='Original Only']
original_ter=landscapes.loc[landscapes['Landscape Category']=='Original and Ter']
reduced=landscapes.loc[landscapes['Landscape Category']=='Reduced Only']
reduced_ter=landscapes.loc[landscapes['Landscape Category']=='Reduced and Ter']

fig, (ax1,ax2,ax3,ax4)=plt.subplots(1,4,sharey=True)
axs=[ax1,ax2,ax3,ax4]
landscape_set=[original, original_ter,reduced, reduced_ter]
names=['Original library','Original library\nand terpolymers','Reduced\noriginal library','Reduced\noriginal library\nand terpolymers']
for i in range(4):
    ax=axs[i]
    landscape=landscape_set[i]
    ax.scatter(landscape['Landscape Number'],landscape['Binary Polymers Score'],c=colors[0],marker=markers[0],label='Binary copolymers')
    ax.scatter(landscape['Landscape Number'],landscape['Terpolymer Test Score'],c=colors[1],marker=markers[1],label='Terpolymers')
    ax.scatter(landscape['Landscape Number'],landscape['Combined Test Score'],c=colors[2],marker=markers[2],label='Combined polymers')
    ax.set_title(names[i],size=11)
    ax.set_xticks([0,1,2,3,4],[1,2,3,4,5])
    ax.tick_params(axis='both',labelsize=11)
ax1.set_ylabel('Test set accuracy',size=11)
ax2.set_xlabel('Number of training sets in landscape',size=11)
fig.set_size_inches(8, 2.5)
plt.savefig('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/Landscape_study.pdf',dpi=300, bbox_inches='tight')
#####make heatmap plots for the increasing data accuracy evaluations
colors=['#debd29','#8bbd31','#527b08','#295208']
cmap1 = LinearSegmentedColormap.from_list("mycmap_kat1",colors)
ogbinary=np.load('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/og_binarytest_increasingdata_eval-rows_og_library_columns_ter.npy')
print(ogbinary)
fig, ax = plt.subplots()
im, cbar = heatmap(ogbinary,[1,2,3,4,5], [1,2,3,4,5],cmap=cmap1,ax=ax,vmin=0.62, vmax=0.96,cbarlabel="Test set accuracy")
for i in range(5):
    for j in range(5):
        text=ax.text(j,i,round(ogbinary[i,j],2),ha="center",va="center",color="black")
fig.tight_layout()
ax.set_xlabel('Number of Terpolymer Training Sets',fontsize=11)
ax.set_ylabel('Number of Binary Copolymer\nTraining Sets',fontsize=11)
ax.set_title(' Test Set',size=11)
fig.set_size_inches(2.5, 2.5)
# plt.savefig('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/ogbinarytest_sized_fixedscale_with_data.pdf',dpi=300, bbox_inches='tight')
#######similarity between the representations
og_only=np.load('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/original_only_similarities.npy')
og_ter=np.load('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/original_ter_similarities.npy')
og_only=np.reshape(og_only,np.size(og_only))
og_ter=np.reshape(og_ter,np.size(og_ter))

fig,(ax1,ax2)=plt.subplots(1,2,sharey=True,sharex=True)
ax2.axvspan(np.mean(og_ter)-np.std(og_ter),np.mean(og_ter)+np.std(og_ter),color=colors[1],alpha=0.2)
ax1.axvspan(np.mean(og_only)-np.std(og_only),np.mean(og_only)+np.std(og_only),color=colors[1],alpha=0.2)
ax1.hist(og_only,color=colors[3],bins=20)
ax2.hist(og_ter,color=colors[2],bins=20)
ax1.set_ylabel('Number of Embedding Pairs',fontsize=11)
ax1.set_xlabel('Cosine Similarity of\nEmbedding Pair',fontsize=11)
ax2.set_xlabel('Cosine Similarity of\nEmbedding Pair',fontsize=11)
ax1.set_title('Only original library polymers')
ax2.set_title('Original Library Polymers \nand Copolymers')
ax1.axvline(np.mean(og_only), color=colors[1], linestyle='dashed', linewidth=1)
# ax1.axvspan(np.mean(og_only)+np.std(og_only),np.mean(og_only)-np.std(og_only),color=colors[1],alpha=0.2)
ax2.axvline(np.mean(og_ter), color=colors[1], linestyle='dashed', linewidth=1)
print('original only',np.mean(og_only),np.median(og_only),np.std(og_only))
print('original and ter',np.mean(og_ter),np.median(og_ter),np.std(og_ter))
# ax2.axvspan(np.mean(og_ter)-np.std(og_ter),np.mean(og_ter)+np.std(og_ter),color=colors[1],alpha=0.2)
plt.savefig('C:/Users/ChemeGrad2020/Documents/Research/Ter-Polymer_Library_Biodeg/Manuscript/similarity_dist.png',dpi=300, bbox_inches='tight')