function [targetVals, inputVals] = targetVals_VS2()
% get input and target values

inputVals.LVm = 120;
inputVals.RVm = 66*2/3; %
inputVals.HR = 180;
inputVals.TBV = 5500;

% those are used to penalize the model during exercise optimization
targetVals.SBP = 180;
targetVals.DBP = 80;
targetVals.EF=0.75;

targetVals.RAPmean = 5; %
targetVals.PPAs = 25;  %
targetVals.PCWP = 15; %
targetVals.CVP = 7;   %
targetVals.CO = 20 ; %
targetVals.sys_dias_time= 0.9 ;

%Information on valves, set to 1 for healthy individuals at rest
%Not used in this form of the model (weight is set to 0)
targetVals.MR = 1; %3
targetVals.MS = 1; %1
targetVals.AI = 1; %1
targetVals.AS = 1; %2 %match bad patient and good and compare parameters
% targetVals.TR = 1; %1
targetVals.TR = 1; %3
targetVals.TS = 1; %1
targetVals.PI = 1; %1
targetVals.PS = 1;%1