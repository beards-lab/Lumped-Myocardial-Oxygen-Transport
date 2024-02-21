# Lumped-Myocardial-Oxygen-Transport
Computer codes associated with Sturgess et al. "Integrated modeling and simulation of recruitment of myocardial perfusion and oxygen delivery in exercise".

# General Overview
This repository contains computer codes associated with Sturgess et al. "Integrated modeling and simulation of recruitment of myocardial perfusion and oxygen delivery in exercise". It contains a whole-body cardiovascular model, a myocardial perfusion model, and an oxygen transport model. The workflow is designed such that results from the whole-body cardiovascular model are used in the myocardial perfusion model and results from the myocardial perfusion model are used in the oxygen transport model. Each computational model is designed to run independently using the required inputs. For user simplicity, required inputs for each model and model results are included in SimulationResults/AllResults.mat 

# Whole Body Cardiovascular Model
The whole body cardiovascular model is run from CVS_TriSeg_Driver.m. Running this code generates pressure waveforms and PV loops for each exercise level using the plotFit.m function. The final section of the code generates the supplementary figures (S1) which show trends of blood pressure, cardiac output, ejection fraction, ventricular volumes, and the ratio of systole to diastole for increasing heart rate. Resting state conditions can be modified by changing values of m and exercise conditions will be modified by changing values of m or a. 
Target values for the resting state are listed in targetVals_Rest and for the maximal exercise state in targetVals_MaxE. Changing values in targetVals_Rest will impact the simulation results as they are used to predict certain parameters and initial conditions.

# Myocardial Perfusion Model
The myocardial perfusion model is run from the MyocardialPerfusionDriver.m. To run this code, you must first run the CVS_TriSegDriver or load SimulationResults/AllResults.mat. This driver only runs for a single exercise level at a given time, so the desired exercise level should be specified in the load data section. Figures showing waveforms of penetrating arterial and venous flow, transmural arterial flow, transmural venous flow, penetrating arterial and venous volumes, transmural arterial volumes, transmural venous volumes, endocardial pressures, and epicardial pressure are generated. 

# Oxygen Transport Model
The oxygen transport model is run from OxygenSimulationDriver.m. To run this code, you must first run the CVS_TriSegDriver or load SimulationResults/ AllResults.mat. This driver only runs for a single exercise level and MVO2 assumption at a given time, so the desired exercise level and MVO2 assumption should be specified in the first two sections.

# PlotMyoModelOxygenResults.m
This code creates a plot of myocardial resistances during different exercise levels and plots related to oxygen saturation and PO2 at the different exercise levels. This should be run after the completion of all 3 models or after opening AllResults.mat



