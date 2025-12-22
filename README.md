Repository for supplementary data and Machine Learning for Structure-Property Relationships for Biodegradability in Copolyesters

Available Supplementary Data:

OD/OD0 plots - found in Biodegradability-OD-Plots

GPC/SEC Traces -found in GPC when available. Subfolders are indicative of solvent system and GPC used for analysis as described in the methods of the published work.

H1 NMR Files analyzed for composition - found in NMR- available as mnova files.

Description of Code Files for Machine Learning
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

  Modeling_Revisions.ipynb
  - used for published modeling results.
  -uses optimized models from original library optimization and reduced library optimization for augmenation with terpolymer data
  -generates splits for the terpolymer library
  -carries out polymer similarity embedding calculations for each landscape used in model training
  -performs retraining of optimized models with increasing amounts of terpolymer data
  -performs updating of optimized models with increasing amounts of terpolymer data
  -extracts feature analysis from optimized models
