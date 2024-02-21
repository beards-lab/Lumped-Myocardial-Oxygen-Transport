%% Set Exercise Level
ExerciseLevel=0;  

if ExerciseLevel==0
    OxygenInput=OxygenInput_Rest;
    RPP_fact=1;
elseif ExerciseLevel==1
    OxygenInput=OxygenInput_MildE1;
    RPP_fact=1.71;
elseif ExerciseLevel==2  
    OxygenInput=OxygenInput_MildE2;
    RPP_fact=2.57;
elseif ExerciseLevel==3
    OxygenInput=OxygenInput_MildE3;
    RPP_fact=3.43;
elseif ExerciseLevel==4  
    OxygenInput=OxygenInput_MaxE;
    RPP_fact=4.30;
end

%% Specify MVO2 hypothesis
% 0= constant transmural MVO2
% 1= transmurally varying MVO2
MVO2_Hyp=0;
%%
C_PA = 0.2/22.4;  %  Molar

i=find(t>tper,1, 'first');

% initial masses
O2x0(1) = C_PA;
O2x0(2) = C_PA;
O2x0(3) = 0.2*C_PA;
O2x0(4) = C_PA;
O2x0(5) = 0.3*C_PA;
O2x0(6) = C_PA;
O2x0(7) = 0.4*C_PA;
O2x0(8) = 0.2*C_PA;
O2x0(9) = 0.25*C_PA;
O2x0(10) = 0.25*C_PA;
O2x0(11) = 0.25*C_PA;


% simulation
 T = t; % this is the time vector from the flow simulation
[O2t,O2x] = ode15s(@dCdT_oxygen,[0 600],O2x0,[],T,RPP_fact,MVO2_Hyp,OxygenInput);
%%
% state variables
O2_output.m_pa = O2x(:,1);
O2_output.m_11 = O2x(:,2);
O2_output.m_21 = O2x(:,3);
O2_output.m_12 = O2x(:,4);
O2_output.m_22 = O2x(:,5);
O2_output.m_13 = O2x(:,6);
O2_output.m_23 = O2x(:,7);
O2_output.m_pv = O2x(:,8);
O2_output.c_t1 = O2x(:,9);
O2_output.c_t2 = O2x(:,10);
O2_output.c_t3 = O2x(:,11);

% volumes

O2_output.v_pa = interp1(T,V_PA,mod(O2t,T(end)));
O2_output.v_11 = interp1(T,V_11,mod(O2t,T(end)));
O2_output.v_12 = interp1(T,V_12,mod(O2t,T(end)));
O2_output.v_21 = interp1(T,V_21,mod(O2t,T(end)));
O2_output.v_22 = interp1(T,V_22,mod(O2t,T(end)));
O2_output.v_13 = interp1(T,V_13,mod(O2t,T(end)));
O2_output.v_23 = interp1(T,V_23,mod(O2t,T(end)));
O2_output.v_pv = interp1(T,V_PV,mod(O2t,T(end))) + 0.1; % add offset to v_pv to keep it stays positive

% concentrations
O2_output.c_pa = O2_output.m_pa./O2_output.v_pa;
O2_output.c_11 = O2_output.m_11./O2_output.v_11;
O2_output.c_21 = O2_output.m_21./O2_output.v_21;
O2_output.c_12 = O2_output.m_12./O2_output.v_12;
O2_output.c_22 = O2_output.m_22./O2_output.v_22;
O2_output.c_13 = O2_output.m_13./O2_output.v_13;
O2_output.c_23 = O2_output.m_23./O2_output.v_23;
O2_output.c_pv = O2_output.m_pv./O2_output.v_pv;

%flows
O2_output.q_pa = interp1(T,Q_PA,mod(O2t,T(end)));
O2_output.q_11 = interp1(T,Q_11,mod(O2t,T(end)));
O2_output.q_m1 = interp1(T,Q_m1,mod(O2t,T(end)));
O2_output.q_12 = interp1(T,Q_12,mod(O2t,T(end)));
O2_output.q_21 = interp1(T,Q_21,mod(O2t,T(end)));
O2_output.q_m2 = interp1(T,Q_m2,mod(O2t,T(end)));
O2_output.q_22 = interp1(T,Q_22,mod(O2t,T(end)));
O2_output.q_13 = interp1(T,Q_13,mod(O2t,T(end)));
O2_output.q_m3 = interp1(T,Q_m3,mod(O2t,T(end)));
O2_output.q_23 = interp1(T,Q_23,mod(O2t,T(end)));
O2_output.q_pv = interp1(T,Q_PV,mod(O2t,T(end)));

%% Calculate oxygen contributions per layer

%epi
t_idx_start=580*tper;
t_idx_end=600*tper;



q_m1_interp=interp1(O2t,O2_output.q_m1,t_idx_start:0.0001:t_idx_end);
O2_output.C_forward_epi=mean(q_m1_interp.*interp1(O2t,O2_output.c_11,t_idx_start:0.0001:t_idx_end).*(q_m1_interp>0));
q_21_interp=interp1(O2t,O2_output.q_21,t_idx_start:0.0001:t_idx_end);
O2_output.C_Back_epi=mean(-q_21_interp.*interp1(O2t,O2_output.c_pv,t_idx_start:0.0001:t_idx_end).*(q_21_interp<0));
O2_output.C_ratio_epi=O2_output.C_Back_epi/O2_output.C_forward_epi;

%mid
q_m2_interp=interp1(O2t,O2_output.q_m2,t_idx_start:0.0001:t_idx_end);
O2_output.C_forward_mid=mean(q_m2_interp.*interp1(O2t,O2_output.c_12,t_idx_start:0.0001:t_idx_end).*(q_m2_interp>0));
q_22_interp=interp1(O2t,O2_output.q_22,t_idx_start:0.0001:t_idx_end);
O2_output.C_Back_mid=mean(-q_22_interp.*interp1(O2t,O2_output.c_pv,t_idx_start:0.0001:t_idx_end).*(q_22_interp<0));
O2_output.C_ratio_mid=O2_output.C_Back_mid/O2_output.C_forward_mid;

%endo
q_m3_interp=interp1(O2t,O2_output.q_m3,t_idx_start:0.0001:t_idx_end);
O2_output.C_forward_endo=mean(q_m3_interp.*interp1(O2t,O2_output.c_13,t_idx_start:0.0001:t_idx_end).*(q_m3_interp>0));
q_23_interp=interp1(O2t,O2_output.q_23,t_idx_start:0.0001:t_idx_end);
O2_output.C_Back_endo=mean(-q_23_interp.*interp1(O2t,O2_output.c_pv,t_idx_start:0.0001:t_idx_end).*(q_23_interp<0));
O2_output.C_ratio_endo=O2_output.C_Back_endo/O2_output.C_forward_endo;

%%
O2t_idx_5per = O2t>570*tper & O2t<590*tper;

O2_output.c_21_max=max(O2_output.c_21(O2t_idx_5per));
O2_output.c_21_min=min(O2_output.c_21(O2t_idx_5per));
O2_output.c_21_avg=mean(interp1(O2t,O2_output.c_21,t_idx_start:0.0001:t_idx_end));

O2_output.c_22_max=max(O2_output.c_22(O2t_idx_5per));
O2_output.c_22_min=min(O2_output.c_22(O2t_idx_5per));
O2_output.c_22_avg=mean(interp1(O2t,O2_output.c_22,t_idx_start:0.0001:t_idx_end));

O2_output.c_23_max=max(O2_output.c_23(O2t_idx_5per));
O2_output.c_23_min=min(O2_output.c_23(O2t_idx_5per));
O2_output.c_23_avg=mean(interp1(O2t,O2_output.c_23,t_idx_start:0.0001:t_idx_end));


%% Save O2 output
if ExerciseLevel==0
    if MVO2_Hyp==0
        O2_Rest=O2_output;
    elseif MVO2_Hyp==1
        O2_Rest_Grad=O2_output;
    end
elseif ExerciseLevel==1
    if MVO2_Hyp==0
        O2_MildE1=O2_output;
    elseif MVO2_Hyp==1
        O2_MildE1_Grad=O2_output;
    end
elseif ExerciseLevel==2  
    if MVO2_Hyp==0
        O2_MildE2=O2_output;
    elseif MVO2_Hyp==1
        O2_MildE2_Grad=O2_output;
    end
elseif ExerciseLevel==3
    if MVO2_Hyp==0
        O2_MildE3=O2_output;
    elseif MVO2_Hyp==1
        O2_MildE3_Grad=O2_output;
    end
elseif ExerciseLevel==4  
    if MVO2_Hyp==0
        O2_MaxE=O2_output;
    elseif MVO2_Hyp==1
        O2_MaxE_Grad=O2_output;
    end
end
%%
t_idx_start=580*tper;
t_idx_end=600*tper;
t_start=580*tper;

figure(106); hold on;
plot(O2t-t_start,O2_output.q_m3.*O2_output.c_13.*(O2_output.q_m3>0), 'DisplayName', 'Arterial Oxygenation Endo', 'LineWidth',2)
plot(O2t-t_start,-O2_output.q_23.*O2_output.c_pv.*(O2_output.q_23<0), 'DisplayName', 'Venous Oxygenation Endo','LineWidth',2)
axis([0 2*tper 0 0.003]); grid
set(gca,'fontsize',15); box on;
title('Origin of Oxygen in SubEndocardium')
ylabel('Oxygen Content * Blood Flow (Molar*mL/s)','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);


 %% Plot of oxygen concentration in distal compartments as function of time
t_start=580*tper;

figure(); clf; axes('position',[0.15 0.15 0.75 0.75]);  hold on;
plot(O2t-t_start,O2_output.c_21,'color','b','linewidth',1.5); 
plot(O2t-t_start,O2_output.c_22,'color','g','linewidth',1.5); 
plot(O2t-t_start,O2_output.c_23,'color','r','linewidth',1.5); 
l = legend('epi','mid','endo'); 
set(l,'fontsize',12,'location','northeast');
set(gca,'fontsize',14); box on;
ylabel('Oxygen Content (Molar)','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);
axis([0 20*tper 0 0.005]); grid
