function [f,O] = dXdT_myocardium_new(t,X,Pdata,m,f_m,flg)

%CONTAINS P_im

% INPUT
P_in = interp1(Pdata(1,:),Pdata(2,:),t);
P_LV = interp1(Pdata(1,:),Pdata(3,:),t);
dPdT = interp1(Pdata(1,:),Pdata(4,:),t);

SIP  = interp1(Pdata(1,:),Pdata(5,:),t);
dSIPdt  = interp1(Pdata(1,:),Pdata(6,:),t);

P_im1 = 0.167*P_LV +SIP;
P_im2 = 0.500*P_LV + SIP;
P_im3 = 0.833*P_LV + SIP;
dPim1_dt = 0.167*dPdT +dSIPdt;
dPim2_dt = 0.500*dPdT +dSIPdt;
dPim3_dt = 0.833*dPdT +dSIPdt;

P_RA=interp1(Pdata(1,:),Pdata(7,:),t);

% PARAMETERS

fact1 = 1*f_m(1);
fact2 = 1*f_m(2);
fact3 = 1*f_m(3);
fact4 = 1*f_m(4);




C_PA = 0.0013*fact4/3;  % mL / mmHg
L_PA = 0.2; % mmHg∙s^2/mL
R_PA = 2.5; % mmHg / (mL / sec)
R_PV = 1.5; % mmHg / (mL / sec)
C_PV = 0.0254/3; % mL / mmHg
R0m = m(5)*44; % mmHg / (mL / sec)
R01 = 1.2*R0m;
R02 = 0.5*R0m;
V01 = m(6)*2.5/3; % mL
V02 = m(7)*8.0/3; % mL
C1 = m(8)*0.013/9; % mL / mmHg
C2 = m(9)*0.254/9; % mL / mmHg
gamma = 0.75; 


%if using ratios from mynard and smolich
cf1 = 0.80*m(1); % epi/endo compliance factor
cf2 = 0.90*m(3); % epi/endo compliance factor
rf1 = 0.80*m(2); % epi/endo resistance factor
rf2 = 0.90*m(4); % mid/endo resistance factor




% STATE VARIABLES
P_PA = X(1); % penetrating artery pressure
Q_PA = X(2); % inlet flow penetrating artery
P11  = X(3); 
P21  = X(4);
P12  = X(5); 
P22  = X(6);
P13  = X(7); 
P23  = X(8);
P_PV = X(9); % penetrating vein pressure

% CALCULATIONS 
V11 = cf1*((P11 - P_im1)*(fact1)*C1+V01);
V21 = cf1*((P21 - P_im1)*C2+V02);
R11 = rf1*R01*(V01/V11)^2;
R21 = rf1*R02*(V02/V21)^2;

Rm1 = R0m*(gamma*R11/R01 + (1-gamma)*R21/R02);
Q11 = (P_PA - P11)/R11;
Qm1 = (P11 - P21)/Rm1;
Q21 = (P21 - P_PV)/R21;
V12 = cf2*((P12 - P_im2)*(fact2)*C1+V01);
V22 = cf2*((P22 - P_im2)*C2+V02);
R12 = rf2*R01*(V01/V12)^2;
R22 = rf2*R02*(V02/V22)^2;
Rm2 = R0m*(gamma*R12/R01 + (1-gamma)*R22/R02);
Q12 = (P_PA - P12)/R12;
Qm2 = (P12 - P22)/Rm2;
Q22 = (P22 - P_PV)/R22;

V13 = (P13 - P_im3)*(fact3)*C1+V01;
V23 = (P23 - P_im3)*C2+V02;
R13 = R01*(V01/V13)^2;
R23 = R02*(V02/V23)^2;
Rm3 = R0m*(gamma*R13/R01 + (1-gamma)*R23/R02);
Q13 = (P_PA - P13)/R13;
Qm3 = (P13 - P23)/Rm3;
Q23 = (P23 - P_PV)/R23;

Q_ima = Q11 + Q12 + Q13;
Q_imv = Q21 + Q22 + Q23;
Q_out = (P_PV - P_RA)/R_PV;

V_PA=(P_PA)*C_PA;
V_PV=(P_PV)*C_PV;

f(1,:) = (Q_PA - Q_ima)/C_PA; % P_PA
f(2,:) = (P_in - P_PA - Q_PA*R_PA)/L_PA; % Q_PA
f(3,:) = (Q11-Qm1)/((fact1)*cf1*C1) + dPim1_dt; % P11
f(4,:) = (Qm1-Q21)/(cf1*C2) + dPim1_dt; % P21
f(5,:) = (Q12-Qm2)/((fact2)*cf2*C1) + dPim2_dt; % P12
f(6,:) = (Qm2-Q22)/(cf2*C2) + dPim2_dt; % P22
f(7,:) = (Q13-Qm3)/((fact3)*C1) + dPim3_dt; % P13
f(8,:) = (Qm3-Q23)/(C2) + dPim3_dt; % P23
f(9,:) = (Q_imv - Q_out)/C_PV; % P_PV

if flag == 0 
    O = []; 
else 
    O= [Q11; Q12; Q13;
    Qm1; Qm2; Qm3;
    Q21; Q22; Q23;
    V11; V12; V13;
    V21; V22; V23;
    P_im1; P_im2; P_im3;
    R11; R12; R13;
    R21; R22; R23;
    Rm1; Rm2; Rm3;
    V_PA; V_PV;
    Q_out];

end 
