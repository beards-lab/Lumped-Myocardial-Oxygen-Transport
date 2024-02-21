function f = dCdT_oxygen(t2,x,T,RPP_fact,MVO2_Hyp,Input)

% if t2>T(end)
%     t=t2-T;
% else 
%     t=t2;
% end
t = mod(t2,T(end));

% volumes
v_pa = interp1(T,Input.V_PA,t);
v_11 = interp1(T,Input.V_11,t);
v_12 = interp1(T,Input.V_12,t);
v_21 = interp1(T,Input.V_21,t);
v_22 = interp1(T,Input.V_22,t);
v_13 = interp1(T,Input.V_13,t);
v_23 = interp1(T,Input.V_23,t); 
v_pv = interp1(T,Input.V_PV,t) + 0.1; % add offset to v_pv to keep it stays positive

%flows
q_pa = interp1(T,Input.Q_PA,t);
q_11 = interp1(T,Input.Q_11,t);
q_m1 = interp1(T,Input.Q_m1,t);
q_12 = interp1(T,Input.Q_12,t);
q_21 = interp1(T,Input.Q_21,t);
q_m2 = interp1(T,Input.Q_m2,t);
q_22 = interp1(T,Input.Q_22,t);
q_13 = interp1(T,Input.Q_13,t);
q_m3 = interp1(T,Input.Q_m3,t);
q_23 = interp1(T,Input.Q_23,t);
q_pv = interp1(T,Input.Q_PV,t);

% masses
m_pa = x(1);
m_11 = x(2);
m_21 = x(3);
m_12 = x(4);
m_22 = x(5);
m_13 = x(6);
m_23 = x(7);
m_pv = x(8);
c_t1 = x(9);
c_t2 = x(10);
c_t3 = x(11);

% concentrations
c_pa = m_pa/v_pa;
c_11 = m_11/v_11;
c_21 = m_21/v_21;
c_12 = m_12/v_12;
c_22 = m_22/v_22;
c_13 = m_13/v_13;
c_23 = m_23/v_23;
c_pv = m_pv/v_pv;

% parameters

CA= 0.2/22.4;  %0.0089 milimoles per mL blood

MVO2=(10.8/22.4/60)*10^-2;  % milimoles of O2 per second per gram of tissue
PS = 100; %permeability surface area product?
v_t=100/3; %33.3 mL per tissue layer


if MVO2_Hyp==0
    % disp("Constant MVO2 between Layers")
    MVO2_endo = RPP_fact*MVO2*v_t;  %mL per sec for whole tissue layer
    MVO2_mid = RPP_fact*MVO2*v_t;
    MVO2_epi = RPP_fact*MVO2*v_t;
else
    R_M=1.5;
%     disp("Varying MVO2 between Layers")
    MVO2_endo = RPP_fact*MVO2*v_t*(2)*(R_M)/(R_M+1);  % 75%extraction-> divide by 3 layers
    MVO2_epi = RPP_fact*MVO2*v_t*(2)*(1)/(R_M+1);  % 75%extraction-> divide by 3 layers
    MVO2_mid = 0.5*(MVO2_endo+MVO2_epi);  % 75%extraction-> divide by 3 layers
end



f(1,:) = q_pa*CA*(q_pa>0) + q_pa*c_pa*(q_pa<0) - q_11*c_pa*(q_11>0) - q_11*c_11*(q_11<0) - ...
                                                 q_12*c_pa*(q_12>0) - q_12*c_12*(q_12<0) - ...
                                                 q_13*c_pa*(q_13>0) - q_13*c_13*(q_13<0); % m_pa
f(2,:) = q_11*c_pa*(q_11>0) + q_11*c_11*(q_11<0) - q_m1*c_11*(q_m1>0) - q_m1*c_21*(q_m1<0); % m_11
f(3,:) = q_m1*c_11*(q_m1>0) + q_m1*c_21*(q_m1<0) - q_21*c_21*(q_21>0) - q_21*c_pv*(q_21<0) ...
            - PS*(c_21 - c_t1); % m_21
f(4,:) = q_12*c_pa*(q_12>0) + q_12*c_12*(q_12<0) - q_m2*c_12*(q_m2>0) - q_m2*c_22*(q_m2<0); % m_12
f(5,:) = q_m2*c_12*(q_m2>0) + q_m2*c_22*(q_m2<0) - q_22*c_22*(q_22>0) - q_22*c_pv*(q_22<0) ...
            - PS*(c_22 - c_t2); % m_22
f(6,:) = q_13*c_pa*(q_13>0) + q_13*c_13*(q_13<0) - q_m3*c_13*(q_m3>0) - q_m3*c_23*(q_m3<0); % m_13
f(7,:) = q_m3*c_13*(q_m3>0) + q_m3*c_23*(q_m3<0) - q_23*c_23*(q_23>0) - q_23*c_pv*(q_23<0) ...
            - PS*(c_23 - c_t3); % m_23
f(8,:) = q_21*c_21*(q_21>0) + q_21*c_pv*(q_21<0) + ...
         q_22*c_22*(q_22>0) + q_22*c_pv*(q_22<0) + ...
         q_23*c_23*(q_23>0) + q_23*c_pv*(q_23<0) - q_pv*c_pv; % m_pv
f(9,:)  = + (PS/v_t)*(c_21 - c_t1) - MVO2_epi/v_t;
f(10,:) = + (PS/v_t)*(c_22 - c_t2) - MVO2_mid/v_t;
f(11,:) = + (PS/v_t)*(c_23 - c_t3) - MVO2_endo/v_t;

