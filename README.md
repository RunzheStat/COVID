# Multi-Objective Model-based Reinforcement Learning for Infectious Disease Control

This repository is the official implementation of the paper "Multi-Objective Model-based Reinforcement Learning for Infectious Disease Control". 

## Requirements
* R version 3.6.3 (2020-02-29)
* Main packages for the proposed estimation procedure
    - abind
    - crayon
    - doFuture
    - doParallel
    - foreach
    - ggplot2
    - ggpubr
    - grid
    - gridExtra
  	- matrixStats
  	- parallel
  	- plyr
  	- randomForest
  	

* Additional packages for experiments
    - gtools
    - dplyr
    - tidyr


## File Overview
1. Files in the main folder: scripts for reproducing results. 

  * `sir_pred.R`: script to conduct estimation and validation of the transition model. 
  * `pareto_eval.R`: script to conduct evaluation of the Pareto optimal policies.
  * `simu_H1N1`: script to conduct evaluation of the Pareto optimal policies for the H1N1 experiment.
  * `simu_valid.R`: script to conduct cross validation of the environment model.
  * `robustness.R`: script to conduct the hyper-parameter sensitivity analysis. 
  * `data.rds`: the data file, consisting of action, observations, and costs of the six cities.
  
  
2. Files in the `/code` folder
    1. the proposed method
        * `SIR.R`: main functions for the generalized dynamics model
        * `RL.R`: main functions for the reinforcement learning algorithms
        * `utility.R`: helper functions for SIR model and RL parts
    2. experiments
        * `simu.R`: main functions for simulation experiments
        * `simu_utility.R`: helper functions for simulation experiments
  
2. Files in the `/SEIR` folder: script to conduct evaluation of the Pareto optimal policies when the DGP is the SEIR model. 

## Reproduce Results

To reproduce ourexperiment results in the paper, open the corresponding script, change the working directory to the main folder which includs this README file, modify the `main_path` in the script, and run commands below. 
