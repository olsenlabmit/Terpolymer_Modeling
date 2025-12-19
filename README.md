Machine Learning Repository for Terpolymer augmentation of a binary copolymer library

Polymer_Similarity.py file contains function definitions from https://github.com/shijiale0609/PolymerEmbedding_EMD for chemical similarity-informed embedding used to generate chemical representations of the polymers

Terpolymer_ML.ipynb:
  -converts bipolymer library splits used in https://doi.org/10.1073/pnas.2220021120 to usable forms for new polymer similarity embedding 
  - generates splits for the reduced original library
  - defines functions for optimizing sckikit-learn Random Forest Classifier Optimization and SDG_Classifier Optimization
  - implements optimization
  - Note: all optimization models are saved as .pkl files

Model_Optimizations.py:
  - carries out polymer embedding calculations and model optimzations using multiprocessing. The same optimizations are can be carried out in Terpolymer_ML.ipynb

Terpolymer_Modeling.ipynb and Terpolymer_Modeling_Updated.ipynb
  -uses optimized models from original library optimization and reduced library optimization for augmenation with terpolymer data
  -generates splits for the terpolymer library
  -carries out polymer similarity embedding calculations for each landscape used in model training
  -performs retraining of optimized models with increasing amounts of terpolymer data
  -performs updating of optimized models with increasing amounts of terpolymer data
