function [params, init] = estimParametersExercise(targets,inputs, modifiers_rest, modifiers_ex_graded)

    %{ 
    Assignment and/or derivation of all model parameters. 
    
    Inputs: 
    targets     - Patient or standard values being fit to 
    inputs      - TBV, wall masses, heart rate
    modifiers   - vector of floats that adjust selected parameters  

    Outputs: 
    params      - vector of parameters used in the model
    init        - initial conditions for ode15s

Originally from https://github.com/beards-lab/CVS_KimRandall by EB Randall, modified by
Filip Jezek

    %} 

    
    
    
    %% Unpack data structure
    
    % we want to use targets and inputs interchangeably. Thus, creating a
    % merged struct. The 'targets' struct is explicitly for cost function.

    inputData = struct();
    ft = fieldnames(targets);
    for i = 1:length(ft)
        inputData.(ft{i}) = targets.(ft{i});
    end
    fi = fieldnames(inputs);
    for i = 1:length(fi)
        inputData.(fi{i}) = inputs.(fi{i});
    end
    
    HR = inputData.HR; % 1/min <-------- Increased for exercise (e.g. 120). Resting: 64 bpm    
    theta=(HR-64)/(180-64);
    
    % Blood pressures (mmHg)
    P_SAs = inputData.SBP; 
    P_SAd = inputData.DBP; 
    PP_sys = P_SAs - P_SAd;
    PP_pul = inputData.PPAs - inputData.PPAd;   %USES REST TARGETS (PPAd)
    
    % Total blood volume (mL) 
    Vtot  = inputData.TBV; 
    
    % Cardiac output (mL s^(-1))
    CO    = inputData.CO / 60*1000;
    
    % End-diastolic and end-systolic pressures (mmHg) and volumes (mL) 
    LVESV = inputData.LVESV;    
    LVEDV = inputData.LVEDV;   RVEDV = inputData.RVEDV;        
    SV=inputData.LVEDV-inputData.LVESV;

    
    %% Volumes
    % snapped at end diastole - maximal ventricles, minimal atria
    
    % Blood volume distribution values - sum total = 1.0 
    d_SA = .15;              d_PA = .05; 
    d_SV = .6;              d_PV = .2;
    Vd = Vtot - LVEDV - RVEDV - inputData.RAVmin - inputData.LAVmin; % distributed volume


    % Total chamber volumes 
    V_SA_0 = d_SA*Vd;      V_PA_0 = d_PA*Vd; 
    V_SV_0 = d_SV*Vd;      V_PV_0 = d_PV*Vd;

    % Unstressed chamber volumes
    V_SA_u = V_SA_0*0.7;   V_PA_u = V_PA_0*0.4; 
    V_SV_u = V_SV_0*0.9;   V_PV_u = V_PV_0*0.9; 
    
    % Stressed chamber volumes
    V_SA_s = V_SA_0 - V_SA_u;  V_PA_s = V_PA_0 - V_PA_u;  
    V_SV_s = V_SV_0 - V_SV_u;  V_PV_s = V_PV_0 - V_PV_u; 
    
%   
                                                
    %% Elastances and compliances 
       
    % Compliances (mL mmHg^(-1))
    C_SA= SV/PP_sys;

    C_SV = V_SV_s/inputData.CVP;

    C_PA=  SV/PP_pul;

    C_PV = V_PV_s/inputData.PCWP; 
    
    %% Resistances 

    % Arteriolar resistances (mmHg s mL^(-1)) 
    MAP = (P_SAs - P_SAd)/3 + P_SAd;
    R_SA = (MAP - inputData.CVP)/CO;  % Systemic arterial resistance
    MPAP = (inputData.PPAs - inputData.PPAd)/3 + inputData.PPAd;
    R_PA = (MPAP - inputData.PCWP)/CO;  % Pulmonary arterial resistance
    
    % Valve resistances (mmHg s mL^(-1)), o for open, c for closed
    % Grade of 1 is used for healthy inividuals, 
    % grade of 2-3 for valve stenosis or reguritation
    if(~isfield(inputData, 'MR') || inputData.MR <= 1.5) 
        R_m_c = Inf; % 0 or 1 -> no regurgitation. 0% RF
    elseif(inputData.MR <= 2.5)
        R_m_c = 1.6; % 2 -> mild. 20% RF
    elseif(inputData.MR <= 3.5)
        R_m_c = .53; % 3 -> moderate. 40% RF
    else
        assert(inputData.MR <= 5 && inputData.MR >= 4, 'Invalid MR Grade Input');
        R_m_c = 0.017; % 4+ -> severe. 60% RF
    end
    
    if(~isfield(inputData, 'MS') || inputData.MS <= 1.5)
        R_m_o = 1.6e-2; % 0 or 1 -> No stenosis. 2.5 mmHg.
    elseif(inputData.MS <= 2.5)
        R_m_o = 2.3e-2; % 2 -> mild. 3.75 mmHg
    elseif(inputData.MS <= 3.5)
        R_m_o = 5.1e-2; % 3 -> moderate. 7.5 mmHg
    else
        assert(inputData.MS <= 5 && inputData.MS >= 4, 'Invalid MS Grade Input');
        R_m_o = 9.4e-2; % 4+ -> severe. 12 mmHg
    end
    
    if(~isfield(inputData, 'AI') || inputData.AI <= 1.5)
        R_a_c = Inf; % 0 or 1 -> No regurgitation. 0% RF
    elseif(inputData.AI <= 2.5)
        R_a_c = 1.8; % 2 -> mild. 20% RF
    elseif(inputData.AI <= 3.5)
        R_a_c = 0.55; % 3 -> moderate. 40% RF
    else 
        assert(inputData.AI <= 5 && inputData.AI >= 4, 'Invalid AI Grade Input');
        R_a_c = 8e-2; % 4+ -> severe. 60% RF
    end

    if(~isfield(inputData, 'AS') || inputData.AS <= 1.5)
        R_a_o = 7.3e-3; % 0 or 1 -> No stenosis. 2.25 mmHg
    elseif(inputData.AS <= 2.5)
        R_a_o = 3.4e-2; % 2 -> mild. 10 mmHg
    elseif(inputData.AS <= 3.5)
        R_a_o = 1.15e-1; % 3 -> moderate. 30 mmHg
    else
        assert(inputData.AS <= 5 && inputData.AS >= 4, 'Invalid AS Grade Input');
        R_a_o = 2.25e-1; % 4+ -> severe. 50 mmHg
    end

    if(~isfield(inputData, 'TR') || inputData.TR <= 1.5)
        R_t_c = Inf; % 0 or 1 -> no regurgitation. 0% RF
    elseif(inputData.TR <= 2.5)
        R_t_c = 1.2e-1; % 2 -> mild. 25% RF
    elseif(inputData.TR <= 3.5)
        R_t_c = 5.2e-2; % 3 -> moderate. 35% RF
    else
        assert(inputData.TR <= 5 && inputData.TR >= 4, 'Invalid TR Grade Input');
        R_t_c = 1.6e-2; % 4+ -> severe. 45% RF
    end

    if(~isfield(inputData, 'TS') || inputData.TS <= 1.5)
        R_t_o = 3e-3; % 0 or 1 -> no stenosis. 0.5 mmHg. 3e-3
    else
        assert(inputData.TS >= 2, 'Invalid TS Grade Input');
        R_t_o = 6e-2; % 2+ -> stenosis. 6 mmHg
    end
    
    if(~isfield(inputData, 'PI') || inputData.PI <= 1.5)
        R_p_c = Inf; % 0 or 1 -> no regurgitation. 0% RF
    elseif(inputData.PI <= 2.5)
        R_p_c = 0.65; % 2 -> mild. 10% RF
    elseif(inputData.PI <= 3.5)
        R_p_c = 1.33e-1; % 3 -> moderate. 30% RF
    else
        assert(inputData.PI <= 5 && inputData.PI >= 4, 'Invalid PI Grade Input');
        R_p_c = 3.3e-2; % 4+ -> severe. 50% RF
    end

    if(~isfield(inputData, 'PS') || inputData.PS <= 1.5)
        R_p_o = 1.3e-3; % 0 or 1 -> no stenosis. 0.67 mmHg
    elseif(inputData.PS <= 2.5)
        R_p_o = 5.85e-2; % 2 -> mild. 22 mmHg
    elseif(inputData.PS <= 3.5)
        R_p_o = 1.7e-1; % 3 -> moderate. 50 mmHg
    else
        assert(inputData.PS <= 5 && inputData.PS >= 4, 'Invalid PS Grade Input');
        R_p_o = 3.82e-1; % 4+ -> severe. 78 mmHg (can't reach that high, so this is doing 76 mmHg). FIXME: breaks when resistance increases more than this (ODE tolerances).
    end
    

    
    %% Heart model parameters 
    
    % Sarcomere length parameters (µm)
    Lsref   = 2; 
    Lsc0    = 1.51; 
    Lse_iso = 0.04; 
    
    % Sarcomere length shortening velocity (µm s^(-1))
    v_max   = .5*7;    
    
    % Passive stress steepness parameter  
    gamma = 7.5; % optimized from ex vivo model 
    
    
    %% Calculate patient-specific reference midwall surface area (Amref) for LV, SEP, and RV
    
    % Ventricular inner chamber radius (cm)
    r_LV_and_SEP = (LVEDV * 3 / (4* pi))^(1/3); 
    r_RV         = (RVEDV * 3 / (4* pi))^(1/3); 
    
    % Ventricle midwall radius (chamber radius (r) + 1/2 wall thickness
    % (H)) (cm)
    H_LW_and_SW = (inputData.Hed_LW + inputData.Hed_SW) / 2; % Not sure which one to use; they're pretty similar
    H_RW         = inputData.Hed_RW / 1; 

    % Midwall radius (cm) 
    r_m_LV_and_SEP = r_LV_and_SEP + H_LW_and_SW/2;  
    r_m_RV         = r_RV + H_RW/2;


    % Outer radius (cm)
    r_o_LV_and_SEP = r_LV_and_SEP + H_LW_and_SW; 
    r_o_RV         = r_RV + H_RW; 
    
    % Midwall reference surface area (cm^2)
    Amref_LV_and_SEP = 4 * pi * (r_m_LV_and_SEP)^2; 
    Am_RV            = 4 * pi * (r_m_RV)^2;
    
    LvSepr = 2/3;
    Amref_LV  = Amref_LV_and_SEP * LvSepr; % Assume LV is 2/3 of LV+SEP 
    Amref_SEP = Amref_LV_and_SEP * (1 - LvSepr); % Assume SEP is 1/3 of LV+SEP
    Amref_RV  = Am_RV;
    
    %% Calculate patient-specific wall volume (Vw) for LV, SEP, and RV 
    
    % Ventricle volume (chamber + midwall) (mL)
    Vw_chamber_LV_and_SEP = 4/3 * pi * r_o_LV_and_SEP^3;  
    Vw_chamber_RV         = 4/3 * pi * r_o_RV^3; 
    
    % Ventricular wall volume (mL)
    Vw_LV_and_SEP = Vw_chamber_LV_and_SEP - LVEDV; 
    Vw_RV         = Vw_chamber_RV - RVEDV;  
    
    % Ventricular midwall volume from data (mL)
    rho_myo = 1.055; % g/mL

    
    Vw_LV  = Vw_LV_and_SEP * LvSepr; % Assume LV is 2/3 of LV+SEP 
    Vw_SEP = Vw_LV_and_SEP * (1 - LvSepr); % Assume SEP is 1/3 of LV+SEP 
    
    %% Approximations for initial displacements and Amref_rv in end-diastole 
    
    % Initialize diastolic displacement values (cm)
    xm_LV_d_0  = -5; 
    xm_SEP_d_0 = 2; 
    xm_RV_d_0  = 6; 
    ym_d_0     = 3; 
    
    x0_d = [xm_LV_d_0; 
        xm_SEP_d_0; 
        xm_RV_d_0;
        ym_d_0; 
        Amref_RV; 
        ]; 
    
    % Inputs for calculating displacements 
    Vw    = [Vw_LV,Vw_SEP,Vw_RV]; 
    Amref = [Amref_LV,Amref_SEP]; 
    
    % Assume end-diastolic sarcomere length 
    SL_d    = 2; %µm 
    
    opts = optimoptions('fsolve','Display','none',...
        'MaxFunctionEvaluations',2e3,'Algorithm','levenberg-marquardt'); 
    [d0,~] = fsolve(@(x) calc_xm_ym(x,Lsref,Vw,Amref,SL_d,LVEDV,0),x0_d,opts); 
    
    
    % Outputs / Diastolic displacements
    xm_LV_d  = d0(1);
    xm_SEP_d = d0(2);
    xm_RV_d  = d0(3);
    ym_d     = d0(4);
    Amref_RV = d0(5); 
    
    % SL lengths (um)
    Lsc_LV_0 = SL_d;
    Lsc_SEP_0 = SL_d;
    Lsc_RV_0 = SL_d;    
    
    init = [xm_LV_d; xm_SEP_d; xm_RV_d; ym_d          % 1-4
        Lsc_LV_0; Lsc_SEP_0; Lsc_RV_0;                % 5-7
        LVEDV; RVEDV; V_SA_s; V_SV_s; V_PA_s; V_PV_s; % 8-13
        inputData.LAVmin; inputData.RAVmin                % 14-15
        ]; 
    
    %% Calculate passive stress parameters (k_pas) for LV and RV in end-diastole  
    
    % Midwall surface area (cm^2)
    Am_LV_d = pi * (xm_LV_d^2  + ym_d^2);
    Am_RV_d = pi * (xm_RV_d^2  + ym_d^2);
    
    % Midwall curvature (cm^(-1))
    Cm_LV_d = 2 * xm_LV_d  / (xm_LV_d^2  + ym_d^2);
    Cm_RV_d = -2 * xm_RV_d  / (xm_RV_d^2  + ym_d^2);
    
    % Midwall ratio (dimensionless) 
    z_LV_d = 3 * Cm_LV_d  * Vw_LV  / (2 * Am_LV_d); 
    z_RV_d = 3 * Cm_RV_d  * Vw_RV  / (2 * Am_RV_d); 
    
    % Instantaneous sarcomere length (µm) in end-diastole
    Ls_LV_d = SL_d; 
    Ls_RV_d = SL_d;  
    
    % Passive stress 
    sigma_pas_LV_d = (Ls_LV_d/Lsc0 - Lsc0/Lsc0)^gamma; 
    sigma_pas_RV_d = (Ls_RV_d/Lsc0 - Lsc0/Lsc0)^gamma; 
    
    % Dimensionless combination function
    Gamma_LV_d = -(2 / 3) * z_LV_d * (1 + (1 / 3) * z_LV_d^2 + (1 / 5) * z_LV_d^4);
    Gamma_RV_d = -(2 / 3) * z_RV_d * (1 + (1 / 3) * z_RV_d^2 + (1 / 5) * z_RV_d^4);
    
    % Passive stress scaling parameters
    k_pas_LV = inputData.PCWP / (Gamma_LV_d * sigma_pas_LV_d);
    k_pas_RV = inputData.CVP / (Gamma_RV_d * sigma_pas_RV_d);
    
    %% Approximations for initial displacements and Amref_rv in end-systole 
    
    % Initialize systolic displacements values (cm)
    xm_LV_s_0  = -5; 
    xm_SEP_s_0 = 2; 
    xm_RV_s_0  = 6; 
    ym_s_0     = 3; 
    
    x0_s = [xm_LV_s_0; 
        xm_SEP_s_0; 
        xm_RV_s_0;
        ym_s_0; 
        Amref_RV; 
        ]; 
    
    Amref = [Amref_LV,Amref_SEP, Amref_RV]; 
    
    opts = optimoptions('fsolve','Display','none',...
        'MaxFunctionEvaluations',2e3,'Algorithm','levenberg-marquardt'); 
    [fnew1,~] = fsolve(@(x) calc_xm_ym(x,Lsref,Vw,Amref,[],LVESV,1),x0_s,opts); 

    
    % Outputs / Systolic displacements
    xm_LV_s = fnew1(1);
    xm_RV_s = fnew1(3);
    ym_s    = fnew1(4);
    
    %% Calculate active stress parameters (k_act) for LV and RV in end-systole 
    
    % Midwall surface area (cm^2)
    Am_LV_s = pi * (xm_LV_s^2  + ym_s^2);
    Am_RV_s = pi * (xm_RV_s^2  + ym_s^2);
    
    % Midwall curvature (cm^(-1))
    Cm_LV_s = 2 * xm_LV_s  / (xm_LV_s^2  + ym_s^2);
    Cm_RV_s = - 2 * xm_RV_s  / (xm_RV_s^2  + ym_s^2);
    
    % Midwall ratio (dimensionless)  
    z_LV_s = 3 * Cm_LV_s  * Vw_LV  / (2 * Am_LV_s); 
    z_RV_s = 3 * Cm_RV_s  * Vw_RV  / (2 * Am_RV_s); 
    
    % Myofiber strain (dimensionless)
    eps_LV_s = 0.5 * log(Am_LV_s  / Amref_LV) - (1/12) * z_LV_s^2  - 0.019 * z_LV_s^4; 
    eps_RV_s = 0.5 * log(Am_RV_s  / Amref_RV) - (1/12) * z_RV_s^2  - 0.019 * z_RV_s^4; 
    
    % Sarcomere length (µm)
    Ls_LV_s  = Lsref * exp(eps_LV_s); 
    Ls_RV_s  = Lsref * exp(eps_RV_s); 
    
    % Activation function 
    Y = .6; % set to 1 in systole 
    
    % Active stress 
    sigma_act_LV_s = Y * (Ls_LV_s/Lsc0  - Lsc0/Lsc0); 
    sigma_act_RV_s = Y * (Ls_RV_s/Lsc0  - Lsc0/Lsc0); 
    
    % Dimensionless combination function 
    Gamma_LV_s = - (2 / 3) * z_LV_s * (1 + (1 / 3) * z_LV_s^2 + (1 / 5) * z_LV_s^4);
    Gamma_RV_s = - (2 / 3) * z_RV_s * (1 + (1 / 3) * z_RV_s^2 + (1 / 5) * z_RV_s^4);
    
    % Active stress scaling parameters 
    k_act_LV = P_SAs / (Gamma_LV_s * sigma_act_LV_s);   
    k_act_RV = inputData.PPAs / (Gamma_RV_s * sigma_act_RV_s);


 

%% Modifying Parameters
 
m = modifiers_rest;
a = modifiers_ex_graded;


%Activation function
% Percentage of cardiac cycle 
params.tau_TS = 0.50*(1-a(8)*theta); % unitless time to maximal systole. Increase leads to longer time for LA in max volume  
params.tau_TR = 0.15*(1-a(8)*theta); % unitless relaxation time (maximal systole to baseline)

% Heart period (s) 
params.HR = HR; % <----------------- Increased in exercise (e.g. 120). Resting: 60 bpm
params.T = 60/HR;

% Compliances (mL mmHg^(-1))

params.C_SA = C_SA*m(1)/(1 + a(1)*theta);  params.C_SV = C_SV*m(15)/(1 + a(2)*theta);  
params.C_PA = C_PA*m(2); params.C_PV = C_PV*m(16);  
    
% Resistances (mmHg s mL^(-1))
params.R_SA  = R_SA*m(3)/(1 + a(3)*theta); % <-------- Scaled down for exercise (1/2). Resting: 1.368
params.R_tSA = 0.08; 
params.R_PA  = R_PA*m(4)/(1 + a(4)*theta); % <-------- Scaled down for exercise (1 / 1.25). Resting: 0.136
params.R_tPA = 0.02; 


params.R_t_o = R_t_o;
params.R_t_c = R_t_c;
params.R_p_o = R_p_o;
params.R_p_c = R_p_c;
params.R_m_o = R_m_o * m(13); % aim for mitral profile
params.R_m_c = R_m_c;   
params.R_a_o = R_a_o;
params.R_a_c = R_a_c;
params.R_RA = 0.040*m(17)/(1 + a(5)*theta);
params.R_LA = 0.040*m(8)/(1 + a(5)*theta); 
    
% Force scaling factors (unitless) 
params.k_pas_LV = k_pas_LV*m(5)/(1 + a(6)*theta);
params.k_pas_RV = k_pas_RV*m(6)/(1 + a(6)*theta);
params.k_act    = k_act_LV*m(7)*((a(7)*theta)+1); % <-------- Scaled up for exercise (1.5). Resting: 2900
    
% Midwall reference surface area (cm^2)
params.Amref_LV  = Amref_LV; 
params.Amref_SEP = Amref_SEP ;  % 1/2 Amref_LV
params.Amref_RV  = Amref_RV*m(12) ; % more than LV
    
% Midwall volume (mL) (Input)
params.Vw_LV  = Vw_LV; %2/3 Total LV Volume (LV FW + SEP)
params.Vw_SEP = Vw_SEP; 
params.Vw_RV  = Vw_RV;

init(13) = init(13)*m(9);

params.RAV0u = inputData.RAVmin*m(10);
params.LAV0u = inputData.LAVmin*m(11);
params.LEa = 2.60 *m(14)*((a(7)*theta)+1);
end 