%% Inputting flow and pressure data
% m = [0.7978    1.0509    0.9604    1.2349    1.3601    0.6105    0.5146    1.5235    1.1013];
m_myo=load("m_MyocardialPerfusion.txt");
%% LOAD DATA 

%specify Exercise Level
%0 for rest; 1 for 90 bpm; 2 for 120 bpm; 3 for 150 bpm; 4 for 180 bpm
ExerciseLevel=0;  
         

% Myo_Targets_Exercise.ENDO_EPI=1.05;
% Myo_Targets_Exercise.ENDO_MID=1.025;

if ExerciseLevel==0
    Myo_Inputs_Hold=Myo_Inputs_Rest;
    HR=64;
elseif ExerciseLevel==1
    Myo_Inputs_Hold=Myo_Inputs_MildE1;
    HR=90;
elseif ExerciseLevel==2  
    Myo_Inputs_Hold=Myo_Inputs_MildE2;
    HR=120;
elseif ExerciseLevel==3
    Myo_Inputs_Hold=Myo_Inputs_MildE3;
    HR=150;
elseif ExerciseLevel==4  
    Myo_Inputs_Hold=Myo_Inputs_MaxE;
    HR=180;
end
    
t=Myo_Inputs_Hold.t;
tstart=t(1);
tend=floor(t(end));
dt=1/1000;
tdata=linspace(tstart,tend,(tend-tstart)/dt);
    
Pao=interp1(t,Myo_Inputs_Hold.P_SA,tdata); 
Plv=interp1(t,Myo_Inputs_Hold.P_LV,tdata);
SIP=interp1(t,Myo_Inputs_Hold.Y, tdata);
SIP=SIP*(0.2*max(Plv)/max(SIP));
Pra=interp1(t,Myo_Inputs_Hold.P_RA,tdata);
tper=60/HR;


dPlvdt(1,:) = (Plv(2)-Plv(1))./dt;
for i = 2:(length(Plv)-1)
  dPlvdt(i,:) = (Plv(i+1)-Plv(i-1))./(2*dt);
end
dPlvdt(length(Plv),:) = (Plv(length(Plv))-Plv(length(Plv)-1))./dt;

dSIPdt(1,:) = (SIP(2)-SIP(1))./dt;
for i = 2:(length(SIP)-1)
  dSIPdt(i,:) = (SIP(i+1)-SIP(i-1))./(2*dt);
end
dSIPdt(length(SIP),:) = (SIP(length(SIP))-SIP(length(SIP)-1))./dt;


Pdata = [tdata; Pao; Plv; dPlvdt'; SIP; dSIPdt'; Pra];


%% Initial Conditions and Simulation

if ExerciseLevel==0
    f=[1,1,1,1];
elseif ExerciseLevel==1
    f=[2.5,3.3,5.0,2.5]; 
elseif ExerciseLevel==2  
    f=[3.7,5.1, 8.7,3.7];
elseif ExerciseLevel==3
    f=[4.9,7.0, 13.0,4.9];
elseif ExerciseLevel==4  
    f=[6.3,9.2,19.2,6.3];
end


Xo_myo = [60 1 50 50 85 85 120 120 5]'; 
opts=odeset('MaxStep',1/200, 'RelTol', 1e-5);

%End time should be set such that it is at the end of a cardiac cycle
%Set to either 30 or 60 seconds or a set number of cardiac cycles 

%Run the ode once to reach steady state 
%and a second time for input data to the oxygen transport model 
[t,X] = ode15s(@dXdT_myocardium,[0.0, 30],Xo_myo,opts,Pdata, m_myo,f,0);  
[t,X] = ode15s(@dXdT_myocardium,[0.0, 10*tper],X(end,:),opts,Pdata, m_myo,f,0);

o = zeros(length(t),30); 
for i = 1:length(t)
    [~,o(i,:)] = dXdT_myocardium(t(i),X(i,:),Pdata,m_myo,f,1); 
end 

Q_PA= X(:,2); 
Q_11  = o(:,1);    
Q_12 = o(:,2);    
Q_13 = o(:,3); 

Q_m1  = o(:,4);  
Q_m2 = o(:,5);    
Q_m3 = o(:,6);

Q_21  = o(:,7);    
Q_22 = o(:,8);    
Q_23 = o(:,9); 

V_11  = o(:,10);    
V_12 = o(:,11);    
V_13 = o(:,12); 

V_21  = o(:,13);    
V_22 = o(:,14);    
V_23 = o(:,15); 

P_im1  = o(:,16);    
P_im2  = o(:,17);    
P_im3  = o(:,18); 

R_11  = o(:,19);    
R_12 = o(:,20);    
R_13 = o(:,21); 

R_21  = o(:,22);    
R_22 = o(:,23);    
R_23 = o(:,24); 

R_m1  = o(:,25);    
R_m2 = o(:,26);    
R_m3 = o(:,27); 

V_PA= o(:,28); 
V_PV= o(:,29);

Q_PV=o(:,30);

PPA  = X(:,1);
P11  = X(:,3); 
P21  = X(:,4);
P12  = X(:,5); 
P22  = X(:,6);
P13  = X(:,7); 
P23  = X(:,8);
PPV  = X(:,9);

volume_epi=V_11+V_21;
volume_mid=V_12+V_22;
volume_endo=V_13+V_23;


t_idx_1per = t>5*tper & t<=10*tper;

Vol_1per=volume_epi(t_idx_1per)+volume_mid(t_idx_1per)+volume_endo(t_idx_1per);

Vol_max=max(Vol_1per);
Vol_min=min(Vol_1per);

t_idx_start=5*tper;
t_idx_end  =10*tper;

Q_epi=mean(interp1(t,Q_11,t_idx_start:0.000001:t_idx_end));
Q_mid=mean(interp1(t,Q_12,t_idx_start:0.000001:t_idx_end));
Q_endo=mean(interp1(t,Q_13,t_idx_start:0.000001:t_idx_end));

QPA_mean=mean(interp1(t,Q_PA,t_idx_start:0.000001:t_idx_end));


avgVolEpi=mean(interp1(t,volume_epi,t_idx_start:0.000001:t_idx_end));
avgVolMid=mean(interp1(t,volume_mid,t_idx_start:0.000001:t_idx_end));
avgVolEndo=mean(interp1(t,volume_endo,t_idx_start:0.000001:t_idx_end));

change_IMBV_Epi=(max(volume_epi(t_idx_1per))-min(volume_epi(t_idx_1per)))/(max(volume_epi(t_idx_1per)));
change_IMBV_Endo=(max(volume_endo(t_idx_1per))-min(volume_endo(t_idx_1per)))/(max(volume_endo(t_idx_1per)));
change_IMBV_Mid=(max(volume_mid(t_idx_1per))-min(volume_mid(t_idx_1per)))/(max(volume_mid(t_idx_1per)));


avgP11=mean(interp1(t,P11,t_idx_start:0.0000001:t_idx_end));
avgP12=mean(interp1(t,P12,t_idx_start:0.0000001:t_idx_end));
avgP13=mean(interp1(t,P13,t_idx_start:0.0000001:t_idx_end));
avgP21=mean(interp1(t,P21,t_idx_start:0.0000001:t_idx_end));
avgP22=mean(interp1(t,P22,t_idx_start:0.0000001:t_idx_end));
avgP23=mean(interp1(t,P23,t_idx_start:0.0000001:t_idx_end));

%Uncomment to display output results
% disp(QPA_mean);
% disp(Q_endo/Q_epi);
% disp(Q_endo/Q_mid);
% 
% disp(avgVolEpi); disp(avgVolMid);disp(avgVolEndo);
% 
% disp((Vol_max-Vol_min)/Vol_max);
%% Save resistance Outputs
t_idx_start=5*tper;
t_idx_end  =10*tper;
t_idx = t>5*tper & t<=10*tper;
idx_first=find(t_idx,1,'first');
idx_last=find(t_idx,1,'last');

    Outputs.QPA_mean=QPA_mean;
    Outputs.EE_ratio=Q_endo/Q_epi;
    Outputs.EM_ratio=Q_endo/Q_mid;
    Outputs.avgVolEpi=avgVolEpi;
    Outputs.avgVolMid=avgVolMid;
    Outputs.avgVolEndo=avgVolEndo;
    Outputs.IMBV_change=(Vol_max-Vol_min)/Vol_max;

    Outputs.R11avg=mean(interp1(t,R_11,t_idx_start:0.0001:t_idx_end));
    Outputs.R12avg=mean(interp1(t,R_12,t_idx_start:0.0001:t_idx_end));
    Outputs.R13avg=mean(interp1(t,R_13,t_idx_start:0.0001:t_idx_end));
    Outputs.Rm1avg=mean(interp1(t,R_m1,t_idx_start:0.0001:t_idx_end));
    Outputs.Rm2avg=mean(interp1(t,R_m2,t_idx_start:0.0001:t_idx_end));
    Outputs.Rm3avg=mean(interp1(t,R_m3,t_idx_start:0.0001:t_idx_end));
    Outputs.R21avg=mean(interp1(t,R_21,t_idx_start:0.0001:t_idx_end));
    Outputs.R22avg=mean(interp1(t,R_22,t_idx_start:0.0001:t_idx_end));
    Outputs.R23avg=mean(interp1(t,R_23,t_idx_start:0.0001:t_idx_end));


    Outputs.R11min=min(R_11(idx_first:idx_last));
    Outputs.R12min=min(R_12(idx_first:idx_last));
    Outputs.R13min=min(R_13(idx_first:idx_last));
    Outputs.Rm1min=min(R_m1(idx_first:idx_last));
    Outputs.Rm2min=min(R_m2(idx_first:idx_last));
    Outputs.Rm3min=min(R_m3(idx_first:idx_last));
    Outputs.R21min=min(R_21(idx_first:idx_last));
    Outputs.R22min=min(R_22(idx_first:idx_last));
    Outputs.R23min=min(R_23(idx_first:idx_last));

    Outputs.R11max=max(R_11(idx_first:idx_last));
    Outputs.R12max=max(R_12(idx_first:idx_last));
    Outputs.R13max=max(R_13(idx_first:idx_last));
    Outputs.Rm1max=max(R_m1(idx_first:idx_last));
    Outputs.Rm2max=max(R_m2(idx_first:idx_last));
    Outputs.Rm3max=max(R_m3(idx_first:idx_last));
    Outputs.R21max=max(R_21(idx_first:idx_last));
    Outputs.R22max=max(R_22(idx_first:idx_last));
    Outputs.R23max=max(R_23(idx_first:idx_last));

%%
%V_PA,V_11,V_12,V_21,V_22,V_13,V_23,V_PV,Q_PA,Q_11,Q_m1,Q_12,Q_21,Q_m2,Q_22,Q_13,Q_m3,Q_23,Q_PV
OxygenInput.V_PA=V_PA;
OxygenInput.V_11=V_11;
OxygenInput.V_12=V_12;
OxygenInput.V_13=V_13;
OxygenInput.V_21=V_21;
OxygenInput.V_22=V_22;
OxygenInput.V_23=V_23;
OxygenInput.V_PV=V_PV;
OxygenInput.Q_PA=Q_PA;
OxygenInput.Q_11=Q_11;
OxygenInput.Q_m1=Q_m1;
OxygenInput.Q_21=Q_21;
OxygenInput.Q_12=Q_12;
OxygenInput.Q_m2=Q_m2;
OxygenInput.Q_22=Q_22;
OxygenInput.Q_13=Q_13;
OxygenInput.Q_m3=Q_m3;
OxygenInput.Q_23=Q_23;
OxygenInput.Q_PV=Q_PV;

if ExerciseLevel==0
    OxygenInput_Rest=OxygenInput;
    MyoModel_Results_Rest=Outputs;
elseif ExerciseLevel==1
    OxygenInput_MildE1=OxygenInput;
    MyoModel_Results_MildE1=Outputs;
elseif ExerciseLevel==2  
    OxygenInput_MildE2=OxygenInput;
    MyoModel_Results_MildE2=Outputs;
elseif ExerciseLevel==3
    OxygenInput_MildE3=OxygenInput;
    MyoModel_Results_MildE3=Outputs;
elseif ExerciseLevel==4  
    OxygenInput_MaxE=OxygenInput;
    MyoModel_Results_MaxE=Outputs;
end


%% Plotting input data
x_start=5*tper;
x_end=10*tper;

figure(); clf; axes('position',[0.15 0.15 0.75 0.75]); hold on;
plot(tdata-x_start,Plv,'k-','linewidth',1.5,'color',0.5*[1 1 1]);
plot(tdata-x_start,Pao,'k-','linewidth',1.5);
plot(tdata-x_start,SIP, 'k-','linewidth',1.5, 'color', 'r')
plot(tdata-x_start,Pra, 'k-','linewidth',1.5, 'color', 'b')
legend('LVP', 'AoP', 'SIP', 'RAP')
set(gca,'Fontsize',15); box on
xlabel('time (sec)','interpreter','latex','fontsize',16);
ylabel('Pressure (mmHg)','interpreter','latex','fontsize',16);
set(gca,'Xlim',[0 2*tper],'Ylim',[0 200]);

%% Plotting Flow
figure(); hold on;
plot(t-x_start,(Q_PA),'color','k','linewidth',2); 
plot(t-x_start,(Q_PV),'color',[0.5,0.5,0.5],'linewidth',2); 
l = legend('arterial','venous'); 
set(l,'fontsize',12,'location','northeast');
set(gca,'fontsize',15); box on;
ylabel('Myocardial flow ','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);
axis([0 2*tper -1 4]); grid
title('Penetrating Artery and Vein Flows')

figure(5); hold on;
plot(t-x_start,Q_13,'color',[180/255,0,0],'linewidth',2); 
plot(t-x_start,Q_12,'color',[0,1,10/25],'linewidth',2); 
plot(t-x_start,Q_11,'color',[0,0,215/255],'linewidth',2); 
l = legend('endo','mid','epi'); 
set(l,'fontsize',12,'location','northeast');
set(gca,'fontsize',15); box on;
title("Arterial Flow")
ylabel('Myocardial flow ','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);
axis([0 2*tper -1 4]); grid

figure(6); hold on;
plot(t-x_start,Q_23,'color',[180/255,0,0],'linewidth',2); 
plot(t-x_start,Q_22,'color',[0,1,10/25],'linewidth',2); 
plot(t-x_start,Q_21,'color',[0,0,215/255],'linewidth',2); 
l = legend('endo','mid','epi'); 
set(l,'fontsize',12,'location','northeast');
set(gca,'fontsize',15); box on;
title('Venous Flow')
ylabel('Myocardial flow ','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);
axis([0 2*tper -1 4]); grid


%% Plotting Pressures 
t_start=x_start;
figure(10); hold on;
plot(t-t_start,PPA,'color',[25/255,210/255,255/255],'linewidth',3, 'DisplayName', "P PA"); 
plot(t-t_start,P11,'color',[130/255,60/255,190/255],'linewidth',3, 'DisplayName', "P 11"); 
plot(t-t_start,P21,'color',[124/255,230/255,4/255],'linewidth',3, 'DisplayName', "P 21"); 
plot(t-t_start,P_im1, 'color', 'k', 'Linestyle', ':','linewidth',3, 'DisplayName', "P IM 1")
plot(t-t_start,PPV,'color',[200/255,100/255,0],'linewidth',3, 'DisplayName', "P PV"); 
set(gca,'fontsize',15); box on;
axis([0 2*tper -50 250]); grid;
title('Epicardial Pressures')
ylabel('Pressure (mmHg)','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);

figure(11); hold on;
plot(t-t_start,PPA,'color',[25/255,210/255,255/255],'linewidth',3, 'DisplayName', "P PA"); 
plot(t-t_start,P13,'color',[130/255,60/255,190/255],'Linestyle', '-','linewidth',3, 'DisplayName', "P 13"); 
plot(t-t_start,P23,'color',[124/255,230/255,4/255],'Linestyle', '-','linewidth',3, 'DisplayName', "P 23"); 
plot(t-t_start,P_im3, 'color', 'k', 'Linestyle', ':','linewidth',3, 'DisplayName', "P IM 3")
plot(t-t_start,PPV,'color',[200/255,100/255,0],'linewidth',3, 'DisplayName', "P PV"); 
axis([0 2*tper -50 250]); grid;
set(gca,'fontsize',15); box on;
title('Endocardial Pressures')
ylabel('Pressure (mmHg)','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);

%% Plotting Volumes
figure(); hold on;
plot(t-x_start,V_13,'color', [180/255,0,0],'linewidth',2, 'DisplayName', 'Endocardial'); 
plot(t-x_start,V_12,'color',[0,1,10/25],'linewidth',2, 'DisplayName', 'Midwall'); 
plot(t-x_start,V_11,'color',[0,0,215/255],'linewidth',2, 'DisplayName', 'Epicardial'); 
axis([0 2*tper 0 3]); grid;
set(gca,'fontsize',15); box on;
title('Proximal Compartment Volumes')
ylabel('Volume (mL)','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);

figure(); hold on;
plot(t-x_start,V_23,'color',[180/255,0,0],'Linestyle', '-','linewidth',2, 'DisplayName', 'Endocardial'); 
plot(t-x_start,V_22,'color',[0,1,10/25],'Linestyle', '-','linewidth',2, 'DisplayName', 'Midwall'); 
plot(t-x_start,V_21,'color',[0,0,215/255],'Linestyle', '-','linewidth',2, 'DisplayName', 'Epicardial'); 
axis([0 2*tper 0 3]); grid;
set(gca,'fontsize',15); box on;
title('Distal Compartment Volumes')
ylabel('Volume (mL)','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);

figure(); hold on;
plot(t-x_start,V_PA,'color','k','Linestyle', '-','linewidth',2, 'DisplayName', 'Penetrating Artery'); 
plot(t-x_start,V_PV,'color',[0.5,0.5,0.5],'Linestyle', '--','linewidth',2, 'DisplayName', 'Penetrating Vein'); 
axis([0 2*tper 0 0.2]); grid;
set(gca,'fontsize',15); box on;
title('Penetrating Artery and Vein Volumes')
ylabel('Volume (mL)','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);


%% Compare Total Myocardial Flow to Patient 3

if ExerciseLevel==0

    %Create 1 wave of model data
    t_idx_2per = t>2*tper & t<=5*tper;
    F_2per=Q_PA(t_idx_2per,1);
    T_2per=t(t_idx_2per,1);
    [pks_A,locs_A]=findpeaks(-F_2per);
    [Min, I]=min(abs(pks_A(2:end)-pks_A(1)));
    
    load("Patient3_UnScaled.mat")
    %Scale Patient Data
    T_Patient3=linspace(0,tper,length(F_Patient));
    Int_Patient=trapz(T_Patient3, F_Patient);
    Flow_Target=QPA_mean*60/HR;
    F_Patient_Scaled3=F_Patient*(Flow_Target/Int_Patient);
    T_Patient3_2=linspace(0,tper*3,3*length(F_Patient));
    
    load("Patient6_UnScaled.mat")
    %Scale Patient Data
    T_Patient6=linspace(0,tper,length(F_Patient));
    Int_Patient=trapz(T_Patient6, F_Patient);
    
    F_Patient_Scaled6=F_Patient*(Flow_Target/Int_Patient);
    T_Patient6_2=linspace(0,3*tper,3*length(F_Patient));
    

    delay=tper-(T_2per(locs_A(1))-T_2per(1));
    figure(30); clf;  hold on;
    plot(T_Patient3_2-delay,[F_Patient_Scaled3,F_Patient_Scaled3,F_Patient_Scaled3],'k', 'linewidth',2);
    plot(T_Patient6_2-delay,[F_Patient_Scaled6,F_Patient_Scaled6,F_Patient_Scaled6],'k','linewidth',2, 'LineStyle',':');
    plot(T_2per(locs_A(1):end)-T_2per(1)-tper,F_2per(locs_A(1):end),'color','r','linewidth',2); 
    l = legend('Patient 3-Scaled','Patient 6-Scaled','Model Flow');
    set(l,'fontsize',12,'location','southeast');
    set(gca,'fontsize',15); box on;
    ylabel('Myocardial flow (mL/sec)','interpreter','latex','fontsize',16);
    xlabel('time (sec)','interpreter','latex','fontsize',16);
    axis([0 2*tper -1 4])
    
    figure(); clf;  hold on;
    plot(tdata-delay,Plv,'k-','linewidth',2,'color',0.5*[1 1 1]);
    plot(tdata-delay,Pao,'k-','linewidth',2);
    plot(tdata-delay,SIP, 'k-','linewidth',2, 'color', 'r')
    plot(tdata-delay,Pra, 'k-','linewidth',2, 'color', 'b')
    legend('LVP', 'AoP', 'SIP', 'RAP')
    set(gca,'Fontsize',15); box on
    xlabel('time (sec)','interpreter','latex','fontsize',16);
    ylabel('Pressure (mmHg)','interpreter','latex','fontsize',16);
    set(gca,'Xlim',[0 2*tper],'Ylim',[0 200]);
end
