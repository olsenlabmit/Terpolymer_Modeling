#EMD calculation from https://github.com/shijiale0609/PolymerEmbedding_EMD.git

#imports --
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

from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split, GridSearchCV, KFold, cross_val_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import pairwise_distances
from scipy.optimize import minimize
from scipy.stats import bernoulli
from scipy.special import expit as sigmoid
from sklearn.datasets import make_moons
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, RationalQuadratic
def EMD_Calculation_pot(
    query_smiles_list=None,
    query_smiles_weight_list=None,
    target_smiles_list=None,
    target_smiles_weight_list=None,
    embedding_function="MorganFingerprint",
    similarity_score_function="Tanimoto",
    radius = 2,
    num_bits = 2048,
):
    # obtain the length of query smiles list and target smiles list
    query_smiles_list_length = len(query_smiles_list)

    target_smiles_list_length = len(target_smiles_list)

    query_smiles_weight_list = list(np.ones([query_smiles_list_length]))

    target_smiles_weight_list = list(np.ones([target_smiles_list_length]))


    # transfer SMILES to fingerprints
    if embedding_function == "RDKFingerprint":

        fpgen = GetRDKitFPGenerator(fpSize=num_bits)

        #fps = [fpgen.GetFingerprint(x) for x in ms]

        query_mol_list = [Chem.MolFromSmiles(x) for x in query_smiles_list]
        query_fingerprint_list = [fpgen.GetFingerprint(x)for x in query_mol_list]
        target_mol_list = [Chem.MolFromSmiles(x) for x in target_smiles_list]
        target_fingerprint_list = [fpgen.GetFingerprint(x) for x in target_mol_list]

    elif embedding_function == "MorganFingerprint":
        query_mol_list = [Chem.MolFromSmiles(x) for x in query_smiles_list]
        query_fingerprint_list = [
            AllChem.GetMorganFingerprintAsBitVect(x, radius, nBits=num_bits)
            for x in query_mol_list
        ]
        target_mol_list = [Chem.MolFromSmiles(x) for x in target_smiles_list]
        target_fingerprint_list = [
            AllChem.GetMorganFingerprintAsBitVect(x, radius, nBits=num_bits)
            for x in target_mol_list
        ]

    elif embedding_function == "MACCSkeys":
        query_mol_list = [Chem.MolFromSmiles(x) for x in query_smiles_list]
        query_fingerprint_list = [MACCSkeys.GenMACCSKeys(x) for x in query_mol_list]
        target_mol_list = [Chem.MolFromSmiles(x) for x in target_smiles_list]
        target_fingerprint_list = [MACCSkeys.GenMACCSKeys(x) for x in target_mol_list]

    else:
        print(
            embedding_function
            + " is not included in the current vision."
            + " Please choose an available embedding function:"
        )
        print("MorganFingerprint, RDKFingerprint, MACCSkeys.")
        return False

    # define the required three sets
    C = np.zeros([query_smiles_list_length,target_smiles_list_length ])

    # use similarity function to calculate d_ij
    if similarity_score_function == "Tanimoto":
        for i in range(0, query_smiles_list_length):
            for j in range(0, target_smiles_list_length):
                C[i,j] = 1 - DataStructs.FingerprintSimilarity(
                    query_fingerprint_list[i], target_fingerprint_list[j],
                    metric=DataStructs.TanimotoSimilarity

                )

    elif similarity_score_function == "Dice":
        for i in range(0, query_smiles_list_length):
            for j in range(0, target_smiles_list_length):
                C[i,j] = 1 - DataStructs.FingerprintSimilarity(
                    query_fingerprint_list[i], target_fingerprint_list[j],
                    metric=DataStructs.DiceSimilarity
                )

    elif similarity_score_function == "Cosine":
        for i in range(0, query_smiles_list_length):
            for j in range(0, target_smiles_list_length):
                C[i,j] = 1 - DataStructs.FingerprintSimilarity(
                    query_fingerprint_list[i], target_fingerprint_list[j],
                    metric=DataStructs.CosineSimilarity
                )

    else:
        print(
            similarity_score_function
            + " is not included in the current vision."
            + " Please choose an available similarity function:"
        )
        print("Tanimoto, Dice, or Cosine")
        return

    #print(C)
    query_smiles_weight_array = np.array(query_smiles_weight_list)/sum(np.array(query_smiles_weight_list))
    target_smiles_weight_array = np.array(target_smiles_weight_list)/sum(np.array(target_smiles_weight_list))
    #print(query_smiles_weight_array, target_smiles_weight_array)
    ot_emd = ot.emd(query_smiles_weight_array, target_smiles_weight_array, C)
    #print(ot_emd)

    W = np.sum(ot_emd * C)
    #print(W)
    return W