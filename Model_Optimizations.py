#author kfransen
#ML study of original and terpolymer librarys

#imports
import numpy as np
import matplotlib.pyplot as plt
import shutil
import sys
import os.path
import json
import seaborn as sns
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.rdFingerprintGenerator import GetRDKitFPGenerator
import ot
import csv
from pickle import dump
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split, GridSearchCV, KFold, cross_val_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import pairwise_distances
from scipy.optimize import minimize
from scipy.stats import bernoulli
from scipy.special import expit as sigmoid
from sklearn.datasets import make_moons
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, RationalQuadratic
from multiprocessing import Pool

#import EMD calculation from Jiale's project
from Polymer_Similarity import EMD_Calculation_pot as EMD
#make some utility functions that end up necessary due to reading csv files or for convenience
def str_to_list(string):
    counter=string.count(',')
    #print(counter)
    list=[]
    for i in range(counter):
        if counter ==1:
            string=string[1:len(string)-1]
        if i==0 and counter>1:
            j=string.index(',')
            list+=[string[2:j-1]]
            string=string[j+1:len(string)-1]
            if string[0]==' ':
                string=string[1:]
        elif i==counter-1:
            j=string.index(',')
            list+=[string[1:j-1]]
            if string[j+1]==' ':
                list+=[string[j+3:len(string)-1]]
            else:
                list+=[string[j+2:len(string)-1]]
        else:
            j=string.index(',')
            list+=[string[1:j-1]]
            if string[0]==' ':
                string=string[j+2:]
            else:
                string=string[j+1:]
        #print(string)

    return list

def concatenate_indices(indices,leave_out):
    #indices is a list of the split arrays, leave out is the index number of the test set
    newlist=[]
    test=None
    for i in range(len(indices)):
        if i==leave_out:
            test=indices[i]
        else:
            newlist+=[indices[i]]
    train=newlist[0]
    for i in range(1,len(newlist)):
        train=np.concatenate((train,newlist[i]),0)
    return [train,test]

def rf_cross_val(C, xs, ys,testx,testy,name):
    #C should be [n_estimators, max_depth,max_leaf_nnodes]
    #xs and ys give a list of the different dataset splices - must be same length, indices in each must match
        #this should be the full training set (will iterate which set is being used as the validation set for cross validation
    scores = []
    for i in range(0, len(xs)):
        X = concatenate_indices(xs, i)
        Y = concatenate_indices(ys, i)
        model = RandomForestClassifier(n_estimators=C[0], max_depth=C[1], min_samples_split=2, min_samples_leaf=1,max_leaf_nodes=C[2])
        model.fit(X[0], Y[0])
        scores += [model.score(X[1], Y[1])]
    score_array = np.array(scores)
    mean = np.mean(score_array)
    std = np.std(score_array)

    #now evaluate against unseen testing data
    X = concatenate_indices(xs, None)
    Y = concatenate_indices(ys, None)
    model = RandomForestClassifier(n_estimators=C[0], max_depth=C[1], min_samples_split=2, min_samples_leaf=1,
                                   max_leaf_nodes=C[2])
    model.fit(X[0], Y[0])
    eval=model.score(testx,testy)
    with open("./240523_Optimization_Models/RF_"+str(name)+"_estim_"+str(C[0])+"_depth_"+str(C[1])+"_leafnodes_"+str(C[2])+".pkl", "wb") as f:
        dump(model, f, protocol=5)
    return [mean, std, eval]

def train_rf(C,xs,ys):
    #C should be [n_estimators, max_depth,max_leaf_nnodes]
    #xs and ys should be list of the training data (no validation set, test set should be left out)
    X = concatenate_indices(xs, i)
    Y = concatenate_indices(ys, i)
    model = RandomForestClassifier(n_estimators=C[0], max_depth=C[1], min_samples_split=2, min_samples_leaf=1,
                                   max_leaf_nodes=C[2])
    model.fit(X[0], Y[0])
    return model

def update_rf(old_model,C,xs,ys):
    #this takes an existing trained rf and fits additional trees given additional training data
    #C should be the number of additional estimators to fit (can be vector)
    #xs and ys should be list of the training data (no validation set, test set should be left out)
    n_estims=len(old_model.estimators_)

    X = concatenate_indices(xs, None)
    Y = concatenate_indices(ys, None)
    model_list=[]
    estimator_amounts=[]
    for i in range(len(C)):
        model = old_model.set_params(n_estimators=C[i]+n_estims,warm_start=True)
        model.fit(X[0], Y[0])
        model_list+=[model]
        estimator_amounts=[len(model.estimators_)]
    return [model_list,estimator_amounts]


def lin_cross_val(C, xs, ys,xtest,ytest,name):
    #C should be [loss, alpha]
    #xs and ys give a list of the different dataset splices - must be same length, indices in each must match
        #this should be the full training set (will iterate which set is being used as the validation set for cross validation
    scores = []
    for i in range(0, len(xs)):
        X = concatenate_indices(xs, i)
        Y = concatenate_indices(ys, i)
        model =SGDClassifier(loss=C[0], alpha=C[1],shuffle=True)
        model.fit(X[0], Y[0])
        scores += [model.score(X[1], Y[1])]
    score_array = np.array(scores)
    mean = np.mean(score_array)
    std = np.std(score_array)
    X = concatenate_indices(xs, None)
    Y = concatenate_indices(ys, None)
    model = SGDClassifier(loss=C[0], alpha=C[1], shuffle=True)
    model.fit(X[0],Y[0])
    eval=model.score(xtest,ytest)
    #save the model
    with open("./240523_Optimization_Models/SGD_"+str(name)+"_loss_"+str(C[0])+"_alpha_"+str(C[1])+".pkl", "wb") as f:
        dump(model, f, protocol=5)
    return [mean, std,eval]

def train_lin(C, xs,ys):
    X = concatenate_indices(xs, i)
    Y = concatenate_indices(ys, i)
    model = SGDClassifier(loss=C[0], alpha=C[1], shuffle=True)
    model.fit(X[0], Y[0])
    return model

def update_lin(old_model, C, xs, ys):
    #xs are additional new data that needs to be fit
    #ys are labels corresponding to xs
    class_biodeg=[0,1] #pretty sure since we use 0 or 1 this is all we've got
    model=old_model.partial_fit(xs,ys,classes=class_biodeg)
    return model

def get_EMD_array(sets,basis,EMD_values=["RDKFingerprint","Tanimoto",5,512]):
    #sets = list of sets (pd dataframes) that need EMD calculations against the basis
    #if basis=None then the concatenation of the sets is used as basis. If not, give polymer basis
    #EMD_values gives [embedding_function, similarity_score_function, radius, num_bits]

    EMD_arrays=[]
    for k in range(len(sets)):
        col = []
        for i in range(len(sets[k].index)):
            row = []
            for j in range(len(basis.index)):
                X = sets[k]['SMILES_List'][i]
                #print(X)
                X = str_to_list(X)
                Y = sets[k]['Ratio_List'][i]
                #print(basis['SMILES_List'])
                #print("basis",basis['SMILES_List'][j])
                Z = basis['SMILES_List'][j]
                #print(Z)
                Z = str_to_list(Z)
                W = basis['Ratio_List'][j]
                # W=str_to_list(W)
                #print(X, Y, Z, W)
                score = EMD(
                    query_smiles_list=X,
                    query_smiles_weight_list=Y,
                    target_smiles_list=Z,
                    target_smiles_weight_list=W,
                    embedding_function=EMD_values[0],
                    similarity_score_function=EMD_values[1],
                    radius=EMD_values[2],
                    num_bits=EMD_values[3],
                )
                row += [score]
            col += [row]
        # save the list of lists as a npy array
        formed = np.array(col)
        EMD_arrays+=[formed]
    return EMD_arrays
#first we need to get the optimization done for the full original library and reduced original library
#load all the data for the original library
print("running the Original EMD calcs")
os0=pd.read_csv('./ML_splits/original_split0.csv',header=0)
os1=pd.read_csv('./ML_splits/original_split1.csv',header=0)
os2=pd.read_csv('./ML_splits/original_split2.csv',header=0)
os3=pd.read_csv('./ML_splits/original_split3.csv',header=0)
os4=pd.read_csv('./ML_splits/original_split4.csv',header=0)
os5=pd.read_csv('./ML_splits/original_split5.csv',header=0)
y_set0=np.load('./ML_splits/Biodeg_fororiginal_split0_vsfull_library.npy')
y_set1=np.load('./ML_splits/Biodeg_fororiginal_split1_vsfull_library.npy')
y_set2=np.load('./ML_splits/Biodeg_fororiginal_split2_vsfull_library.npy')
y_set3=np.load('./ML_splits/Biodeg_fororiginal_split3_vsfull_library.npy')
y_set4=np.load('./ML_splits/Biodeg_fororiginal_split4_vsfull_library.npy')
y_set5=np.load('./ML_splits/Biodeg_fororiginal_split5_vsfull_library.npy')
os_train=[os0,os1,os2,os3,os4] #the training data used for cross validation
#print("concatenated training")
#print(pd.concat(os_train,axis=0).reset_index())
os_train_emd=get_EMD_array(os_train,pd.concat(os_train).reset_index())
os_train_ys=[y_set0,y_set1,y_set2,y_set3,y_set4]
os_test_emd=get_EMD_array([os5],pd.concat(os_train).reset_index())
os_test_ys=y_set5
###### now load and set up the EMD scoring for the reduced library

print("running the reduced EMD calcs")
rl0=pd.read_csv('./Reduced_ML_Splits/reduced_split0.csv',header=0)
rl1=pd.read_csv('./Reduced_ML_Splits/reduced_split1.csv',header=0)
rl2=pd.read_csv('./Reduced_ML_Splits/reduced_split2.csv',header=0)
rl3=pd.read_csv('./Reduced_ML_Splits/reduced_split3.csv',header=0)
rl4=pd.read_csv('./Reduced_ML_Splits/reduced_split4.csv',header=0)
rl5=pd.read_csv('./Reduced_ML_Splits/reduced_split5.csv',header=0)

rl_train=[rl0,rl1,rl2,rl3,rl4] #the training data used for cross validation
rl_train_emd=get_EMD_array(rl_train,pd.concat(rl_train).reset_index())
#load the y values
rl_train_ys=[]
for i in range(5):
    sample=rl_train[i]
    nparray=sample["Biodegradability"].to_numpy()
    rl_train_ys+=[nparray]
rl_test_emd=get_EMD_array([rl5],pd.concat(rl_train).reset_index())
rl_test_ys=rl5["Biodegradability"].to_numpy()


#optimize across different parameters
number_estimators = [2**0, 2**1, 2**2, 2**3, 2**4, 2**5, 2**6, 2**7, 2**8, 2**9, 2**10]
max_depth = [1,2,4,8,16,32,64,128,256]
min_samples_split = [2]
min_samples_leaf = [1]
max_leaf_nodes = [2,4,8,16,32,64,128,256,512]
#make a list of all the combinations of parameters to run
combos_rf=[]
for i in number_estimators:
    for j in max_depth:
        for k in max_leaf_nodes:
            combos_rf+=[[i,j,k]]
#print(combos)

#get the trained values for full library and write them into a file to save parameters, cross val score, and error
# with open('original_library_training_optimization_rf.csv', 'w', newline='') as file:
#         file_writer=csv.writer(file, delimiter=',')
#         file_writer.writerow(['Number of Estimators','Max Depth','Min Samples Split','Min Samples Leaf','Max Leaf Nodes','Mean Cross Val Score','STD','Test Evaluation Score'])
# #get the trained values and write them into a file to save parameters, cross val score, and error
# with open('reduced_library_training_optimization_rf.csv', 'w', newline='') as file:
#         file_writer=csv.writer(file, delimiter=',')
#         file_writer.writerow(['Number of Estimators','Max Depth','Min Samples Split','Min Samples Leaf','Max Leaf Nodes','Mean Cross Val Score','STD','Test Evaluation Score'])


#make this efficient by using multiprocessing modules to maximize CPU usage
def task(arg):
    com=arg
    #execute across original library
    output=rf_cross_val(com,os_train_emd,os_train_ys,os_test_emd[0],os_test_ys,'os')
    with open('original_library_training_optimization_rf.csv', 'a', newline='') as file:
        file_writer=csv.writer(file, delimiter=',')
        file_writer.writerow([com[0],com[1],2,1,com[2],output[0],output[1],output[2]])
    #execute across reduced library
    output = rf_cross_val(com, rl_train_emd, rl_train_ys, rl_test_emd[0], rl_test_ys,"rl")
    with open('reduced_library_training_optimization_rf.csv', 'a', newline='') as file:
        file_writer = csv.writer(file, delimiter=',')
        file_writer.writerow([com[0], com[1], 2, 1, com[2], output[0], output[1], output[2]])
    print('done combo ',str(com),flush=True)


#set up the SGD combinations that we are interested in optimizing over
losses=['hinge','log_loss','perceptron']
alphas=[0.0001,0.001,0.01,0.1,1,10,100]

combos_SDG=[]
for i in losses:
    for j in alphas:
        combos_SDG+=[[i,j]]

#get the trained values for full library and write them into a file to save parameters, cross val score, and error
with open('original_library_training_optimization_SDG.csv', 'w', newline='') as file:
        file_writer=csv.writer(file, delimiter=',')
        file_writer.writerow(['Loss Function','Alpha','Mean Cross Val Score','STD','Test Evaluation Score'])
#get the trained values and write them into a file to save parameters, cross val score, and error
with open('reduced_library_training_optimization_SDG.csv', 'w', newline='') as file:
        file_writer=csv.writer(file, delimiter=',')
        file_writer.writerow(['Loss Function','Alpha','Mean Cross Val Score','STD','Test Evaluation Score'])

def task2(arg):
    com=arg
    # execute across original library
    output = lin_cross_val(com, os_train_emd, os_train_ys, os_test_emd[0], os_test_ys,'os')
    with open('original_library_training_optimization_SDG.csv', 'a', newline='') as file:
        file_writer = csv.writer(file, delimiter=',')
        file_writer.writerow([com[0], com[1], output[0], output[1], output[2]])
    # execute across reduced library
    output = lin_cross_val(com, rl_train_emd, rl_train_ys, rl_test_emd[0], rl_test_ys,'rl')
    with open('reduced_library_training_optimization_SDG.csv', 'a', newline='') as file:
        file_writer = csv.writer(file, delimiter=',')
        file_writer.writerow([com[0], com[1], output[0], output[1], output[2]])
    print('done combo ', str(com), flush=True)
print("made it to multiprocessing start")
if __name__ =='__main__':
    #with Pool() as pool:
     #   pool.map(task,combos_rf)
    with Pool() as pool:
        pool.map(task2,combos_SDG)




