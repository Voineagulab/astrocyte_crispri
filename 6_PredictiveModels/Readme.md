# Repository Contents

This repository contains the data and analysis scripts corresponding to the section on **Benchmarking and improving enhancer-gene regulatory prediction models** in Green *et al.* 2024.

## Data

- **InputData**: Data used to assemble the training data for astrocytes and K562 cells. Details on each file are listed in `FileDescriptions.txt`.
- **TrainingData**: 
  - `TrainingDataframe_Astrocytes.csv`
  - `TrainingDataframe_K562.csv`.

## Scripts

- `GenerateTrainingData_Astrocytes.R`: Generates `TrainingDataframe_Astrocytes.csv`, using data in `/InputData/Astrocytes`.
- `GenerateTrainingData_K562.R`: Generates `TrainingDataframe_K562.csv`, using data in `/InputData/K562`.
- `RunRFmodels.R`: Carries out all the random forest model training and prediction.
- `Functions.R`: Functions sourced in `RunRFmodels.R`.
- `EvaluatePredictionModels.R`: Carries out the model performance evaluations and produces the manuscript figures and supplementary tables.
- `EvaluatePredictionModelFunctions.R`: Functions sourced in `EvaluatePredictionModels.R`.
- `EGrf_AllIntergenicPeaks_Pred.R`: Script that applies EGrf to all intergenic peaks not tested in the CRISPRi screen and predicts their effect on genes within 500kb.

## Predictions

- **ENCODE_Predictions**: Enhancer-gene regulatory prediction data obtained from the ENCODE consortium for Astrocytes (rE2G, ABC models) and K562 cells (rE2G, rE2G-extended, ABC models).
- **RF_Results**: Contains the output of `RunRFmodels.R` for astrocytes and K562 cells. The main results files are provided:
  - `Astrocyte_trainingData_pluspredictions.csv`
  - `K562_trainingData_pluspredictions.csv`. 
  Running the script will also generate additional files; specifically, for each RF model: the hyper parameters, variable importance, and the RF model with cross-validation in `.rds` format.

## Figures and Tables

- **Figure 8 plots**: Generated by `EvaluatePredictionModels.R`. Running the script will also generate supplementary figures and additional plots.