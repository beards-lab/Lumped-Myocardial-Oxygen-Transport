function [targetVals, inputVals] = targetVals_VS1()
% get input and target values
% These input and target parameters are used to predict initial parameter
% values, changing these values will change the simulation results

% male
inputVals.LVm = 120;
inputVals.RVm = 66*2/3; 
inputVals.HR = 64;
inputVals.TBV = 5500;

% those are used to penalize the model during optimizing
targetVals.SBP = 120;
targetVals.DBP = 80;
targetVals.LVEDV = 150;
targetVals.LVESV = 60;
targetVals.Hed_LW = 0.93; % * Posterior LW thickness, end-diastole (cm)
targetVals.Hed_SW = 0.92; % * SW thickness, end-diastole (cm)
targetVals.EAr = 1.7;
targetVals.LAVmax = 70;
targetVals.LAVmin = 30;
targetVals.RVEDV = 165;
targetVals.RVESV = 75;
targetVals.Hed_RW = 0.35; % RW thickness, end-diastole (cm)
targetVals.RAVmax = 70;
targetVals.RAVmin = 30;
targetVals.RAPmax = 6; %*
targetVals.RAPmin = 2; %*
targetVals.PPAs = 22.5;  %
targetVals.PPAd = 11.5;  %
targetVals.PCWP = 9; %*
targetVals.CVP = 4;   %*
targetVals.CO = 5.76 ;
targetVals.sys_dias_time= 0.65;



%Information on valves, set to 1 for healthy individuals 
%Not used in this form of the model (weight is set to 0)
targetVals.MR = 1; %
targetVals.MS = 1; %
targetVals.AI = 1; %
targetVals.AS = 1; %
% targetVals.TR = 1; %
targetVals.TR = 1; %
targetVals.TS = 1; %
targetVals.PI = 1; %
targetVals.PS = 1;%