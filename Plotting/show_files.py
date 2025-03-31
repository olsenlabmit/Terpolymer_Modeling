# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 16:33:11 2024

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
from glob import glob


data_folder = 'E:/Terpolymer Library Data/Biodegradation Images/Last_Batch/'#Analyzed Polymers Keep/'

date= '20231209'


samples = sorted(glob(data_folder + '{}_*'.format(date)))
init=0
for s, sample in enumerate(samples[init:]):
    print(s,sample)
    try:
        try:
            print(np.load(sample + '/biodegradation_classification.npy'))
        except:
            print(np.load(sample + '/classified_by_rate.npy'))
        img=mpl.image.imread(sample + '/plot_wellplate_od_limit_divOD0_cleaned.png')
        plt.imshow(img)
        plt.show()
    except:
        print('no image')