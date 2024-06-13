% A standalone endpoint to run the TriSeg Model 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Setting Run Information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


baseline_HR=64;
Max_HR=180;

Scaled_HR1=90;
Scaled_HR2=120;
Scaled_HR3=150;




%% Run for rest

[targets_rest, inputs]= targetVals_Rest();
inputs.HR=baseline_HR;
[params, init] = estimParametersExercise_graded(targets_rest,inputs);

targets=targets_rest;

printStats = true;
runSimExercise
WholeBody_Results_Rest=o;
Myo_Inputs_Rest.P_SA=P_SA;
Myo_Inputs_Rest.P_LV=P_LV;
Myo_Inputs_Rest.P_RA=P_RA;
Myo_Inputs_Rest.Y=Y;
Myo_Inputs_Rest.t=t;
plotFit


%%
%Run for max exercise
[targets_rest, inputs]= targetVals_Rest();
inputs.HR=Max_HR;
[params, init] = estimParametersExercise_graded(targets_rest,inputs);

% Uncomment the following line to compare max exercise to the max exercise
% targets
%[targets, ~]= targetVals_VS2();

printStats = true;
runSimExercise
WholeBody_Results_MaxE=o;
Myo_Inputs_MaxE.P_SA=P_SA;
Myo_Inputs_MaxE.P_LV=P_LV;
Myo_Inputs_MaxE.P_RA=P_RA;
Myo_Inputs_MaxE.Y=Y;
Myo_Inputs_MaxE.t=t;
plotFit

%%
% Run for mild exercise 1 HR=90 bpm
[targets_rest, inputs]= targetVals_Rest();
inputs.HR=Scaled_HR1;
[params, init] = estimParametersExercise_graded(targets_rest,inputs);

printStats = true;
targets;
runSimExercise
WholeBody_Results_MildE1=o;
Myo_Inputs_MildE1.P_SA=P_SA;
Myo_Inputs_MildE1.P_LV=P_LV;
Myo_Inputs_MildE1.P_RA=P_RA;
Myo_Inputs_MildE1.Y=Y;
Myo_Inputs_MildE1.t=t;
%plotFit

%%
%Run for mild exercise 2 HR=120 bpm
[targets_rest, inputs]= targetVals_Rest();
inputs.HR=Scaled_HR2;
[params, init] = estimParametersExercise_graded(targets_rest,inputs);

printStats = true;
runSimExercise
WholeBody_Results_MildE2=o;
Myo_Inputs_MildE2.P_SA=P_SA;
Myo_Inputs_MildE2.P_LV=P_LV;
Myo_Inputs_MildE2.P_RA=P_RA;
Myo_Inputs_MildE2.Y=Y;
Myo_Inputs_MildE2.t=t;
plotFit

%%
%Run for mild exercise 3 HR=150 bpm
[targets_rest, inputs]= targetVals_Rest();
inputs.HR=Scaled_HR3;
[params, init] = estimParametersExercise_graded(targets_rest,inputs);

printStats = true;
runSimExercise
WholeBody_Results_MildE3=o;
Myo_Inputs_MildE3.P_SA=P_SA;
Myo_Inputs_MildE3.P_LV=P_LV;
Myo_Inputs_MildE3.P_RA=P_RA;
Myo_Inputs_MildE3.Y=Y;
Myo_Inputs_MildE3.t=t;
plotFit

%% Plotting supplemental figures S1 A-F
%To run this section, either run the above code 
% or load WholeBodyResults from SimulationResults
figure(); hold on;
plot([64,90,120,150,180],[WholeBody_Results_Rest.SBP,WholeBody_Results_MildE1.SBP,WholeBody_Results_MildE2.SBP,WholeBody_Results_MildE3.SBP,WholeBody_Results_MaxE.SBP], '-o', 'LineWidth', 3, 'DisplayName','SBP')
plot([64,90,120,150,180],[WholeBody_Results_Rest.DBP,WholeBody_Results_MildE1.DBP,WholeBody_Results_MildE2.DBP,WholeBody_Results_MildE3.DBP,WholeBody_Results_MaxE.DBP], '-o', 'LineWidth', 3, 'DisplayName','DBP')
title('Pressures During Exercise')
xlabel("Heart Rate (BPM)")
ylabel("Pressure (mmHg)")
legend
ylim([0,185])
set(gca,'FontSize',15,'LineWidth',1,'TickDir','both','TickLength',[0.01 0.05]);

figure(); hold on;
plot([64,90,120,150,180],[WholeBody_Results_Rest.EF,WholeBody_Results_MildE1.EF,WholeBody_Results_MildE2.EF,WholeBody_Results_MildE3.EF,WholeBody_Results_MaxE.EF], '-o','LineWidth', 3, 'DisplayName','EF')
title('Ejection Fraction During Exercise')
xlabel("Heart Rate (BPM)")
ylabel("Ejection Fraction (%)")
set(gca,'FontSize',15,'LineWidth',1,'TickDir','both','TickLength',[0.01 0.05]);
ylim([0,0.75])

figure(); hold on;
plot([64,90,120,150,180],[WholeBody_Results_Rest.CO,WholeBody_Results_MildE1.CO,WholeBody_Results_MildE2.CO,WholeBody_Results_MildE3.CO,WholeBody_Results_MaxE.CO], '-o', 'LineWidth', 3 )
set(gca,'FontSize',15,'LineWidth',1,'TickDir','both','TickLength',[0.01 0.05]);
title('Cardiac Output During Exercise')
xlabel("Heart Rate (BPM)")
ylabel("Cardiac Output (L per min)")
ylim([0,21])

figure(); hold on;
plot([64,90,120,150,180],[WholeBody_Results_Rest.sys_dias_time, WholeBody_Results_MildE1.sys_dias_time,WholeBody_Results_MildE2.sys_dias_time,WholeBody_Results_MildE3.sys_dias_time, WholeBody_Results_MaxE.sys_dias_time], '-o','LineWidth', 3)
ylim([0,1])
set(gca,'FontSize',15,'LineWidth',1,'TickDir','both','TickLength',[0.01 0.05]);
title('Systole to Diastole Time During Exercise')
xlabel("Heart Rate (BPM)")
ylabel("Systole To Diastole Time Ratio")
ylim([0,1])

figure(); hold on;
plot([64,90,120,150,180],[WholeBody_Results_Rest.LVEDV,WholeBody_Results_MildE1.LVEDV,WholeBody_Results_MildE2.LVEDV, WholeBody_Results_MildE3.LVEDV, WholeBody_Results_MaxE.LVEDV], '-o', 'LineWidth', 3,'DisplayName','LV EDV')
plot([64,90,120,150,180],[WholeBody_Results_Rest.LVESV,WholeBody_Results_MildE1.LVESV,WholeBody_Results_MildE2.LVESV, WholeBody_Results_MildE3.LVESV, WholeBody_Results_MaxE.LVESV], '-o', 'LineWidth', 3,'DisplayName','LV ESV')
set(gca,'FontSize',15,'LineWidth',1,'TickDir','both','TickLength',[0.01 0.05]);
title('Left Ventricle Volumes')
xlabel("Heart Rate (BPM)")
ylabel("Volume (mL)")
legend
ylim([0,180])



figure(); hold on;
plot([64,90,120,150,180],[WholeBody_Results_Rest.LVEDV-WholeBody_Results_Rest.LVESV, WholeBody_Results_MildE1.LVEDV-WholeBody_Results_MildE1.LVESV, ...
   WholeBody_Results_MildE2.LVEDV-WholeBody_Results_MildE2.LVESV, WholeBody_Results_MildE3.LVEDV-WholeBody_Results_MildE3.LVESV, WholeBody_Results_MaxE.LVEDV- WholeBody_Results_MaxE.LVESV], '-o','LineWidth', 3)
set(gca,'FontSize',15,'LineWidth',1,'TickDir','both','TickLength',[0.01 0.05]);
title('Left Ventricle Stroke Volume')
xlabel("Heart Rate (BPM)")
ylabel("Volume (mL)")
ylim([0,120])


