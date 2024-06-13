% Runs the simulation and gathers all outputs.
% Can be a function, but is a script. whatever.

%% Simulations 
T = params.T;
HR = params.HR;
M = eye(length(init));
M(1,1) = 0;
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0; 
options = odeset('Mass',M,'RelTol',1e-8,'AbsTol',1e-8,'MaxStep',T/100);

% Assume steady state in 20 then run for 30 seconds total
% 30 secons represents a whole number of cardiac cycles for each exercise
% level

[t,y]  = ode15s(@dXdT,[0 20*T],init,options, params);
[t,y]  = ode15s(@dXdT,[0 30],y(end,:),options, params);


% Outputs 
o = zeros(45,length(t));  
for i = 1:length(t) 
  [~,o(:,i)] = dXdT(t(i),y(i,:), params);
end 
   
xm_LV  = y(:,1); 
xm_SEP = y(:,2); 
xm_RV  = y(:,3); 
ym     = y(:,4);     
Lsc_LV  = y(:,5);
Lsc_SEP = y(:,6); 
Lsc_RV  = y(:,7);     
V_LV = y(:,8); 
V_RV = y(:,9); 
V_SA = y(:,10); 
V_SV = y(:,11); 
V_PA = y(:,12); 
V_PV = y(:,13); 
Vtot = sum(y(end,8:13)) ;     
P_LV = o(1,:)'; 
P_SA = o(2,:)'; 
P_SV = o(3,:)'; 
P_RV = o(4,:)'; 
P_PA = o(5,:)'; 
P_PV = o(6,:)';    
Vm_LV  = o(7,:)'; 
Vm_SEP = o(8,:)'; 
Vm_RV  = o(9,:)';     
Am_LV  = o(19,:)'; 
Am_SEP = o(11,:)'; 
Am_RV  = o(12,:)';     
Cm_LV  = o(13,:)'; 
Cm_SEP = o(14,:)';
Cm_RV  = o(15,:)';     
eps_LV  = o(16,:)'; 
eps_SEP = o(17,:)'; 
eps_RV  = o(18,:)'; 
sigma_pas_LV  = o(19,:)';
sigma_pas_SEP = o(20,:)';
sigma_pas_RV  = o(21,:)';    
sigma_act_LV  = o(22,:)';
sigma_act_SEP = o(23,:)';
sigma_act_RV  = o(24,:)';
sigma_LV  = o(25,:)';
sigma_SEP = o(26,:)';
sigma_RV  = o(27,:)';
Q_m = o(28,:)' ;    % Flow across mitral valve (QIN_LV)
Q_a = o(29,:)';     % Flow across aortic valve (QOUT_LV)
Q_t = o(30,:)' ;    % Flow across tricuspid valve (QIN_RV)
Q_p = o(31,:)' ;    % Flow across pulmonary valve (QOUT_RV)
Q_SA = o(32,:)' ; 
Q_PA = o(33,:)' ;    
Tm_LV  = o(34,:)';
Tm_SEP = o(35,:)'; 
Tm_RV  = o(36,:)'; 
Y = o(37,:)'; 
V_RA = o(38,:)';
V_LA = o(39,:)';
P_RA = o(40,:)';
P_LA = o(41,:)';
QIN_RA = o(42,:)';
d_LW = o(43, :)';
d_SW = o(44, :)';
d_RW = o(45, :)';

%% Figures

% FLOW INDICES
% Mitral
Qm_sign = sign(Q_m);
if(Qm_sign(1) <= 0)
    Qm_pos_start = find(Qm_sign == 1, 1);
    Qm_sign = Qm_sign(Qm_pos_start: end);
    Qm_pos_end = find(Qm_sign ~= 1, 1) + Qm_pos_start - 2;
    Qm_neg_start = Qm_pos_end + 1;
    Qm_sign = Qm_sign(Qm_neg_start - Qm_pos_start + 1: end);
    Qm_neg_end = find(Qm_sign == 1, 1) + Qm_neg_start - 2;
else
    Qm_neg_start = find(Qm_sign ~= 1, 1);
    Qm_sign = Qm_sign(Qm_neg_start: end);
    Qm_neg_end = find(Qm_sign == 1, 1) + Qm_neg_start - 2;
    Qm_pos_start = Qm_neg_end + 1;
    Qm_sign = Qm_sign(Qm_pos_start - Qm_neg_start + 1: end);
    Qm_pos_end = find(Qm_sign ~= 1, 1) + Qm_pos_start - 2;
end
Qm_pos = [Qm_pos_start: Qm_pos_end]; % Indices for positive mitral flow
Qm_neg = [Qm_neg_start: Qm_neg_end]; % Indices for negative mitral flow
% Aortic
Qa_sign = sign(Q_a);
if(Qa_sign(1) <= 0)
    Qa_pos_start = find(Qa_sign == 1, 1);
    Qa_sign = Qa_sign(Qa_pos_start: end);
    Qa_pos_end = find(Qa_sign ~= 1, 1) + Qa_pos_start - 2;
    Qa_neg_start = Qa_pos_end + 1;
    Qa_sign = Qa_sign(Qa_neg_start - Qa_pos_start + 1: end);
    Qa_neg_end = find(Qa_sign == 1, 1) + Qa_neg_start - 2;
else
    Qa_neg_start = find(Qa_sign ~= 1, 1);
    Qa_sign = Qa_sign(Qa_neg_start: end);
    Qa_neg_end = find(Qa_sign == 1, 1) + Qa_neg_start - 2;
    Qa_pos_start = Qa_neg_end + 1;
    Qa_sign = Qa_sign(Qa_pos_start - Qa_neg_start + 1: end);
    Qa_pos_end = find(Qa_sign ~= 1, 1) + Qa_pos_start - 2;
end
Qa_pos = [Qa_pos_start: Qa_pos_end]; % Indices for positive aortic flow
Qa_neg = [Qa_neg_start: Qa_neg_end]; % Indices for negative aortic flow
% Tricuspid
Qt_sign = sign(Q_t);
if(Qt_sign(1) <= 0)
    Qt_pos_start = find(Qt_sign == 1, 1);
    Qt_sign = Qt_sign(Qt_pos_start: end);
    Qt_pos_end = find(Qt_sign ~= 1, 1) + Qt_pos_start - 2;
    Qt_neg_start = Qt_pos_end + 1;
    Qt_sign = Qt_sign(Qt_neg_start - Qt_pos_start + 1: end);
    Qt_neg_end = find(Qt_sign == 1, 1) + Qt_neg_start - 2;
else
    Qt_neg_start = find(Qt_sign ~= 1, 1);
    Qt_sign = Qt_sign(Qt_neg_start: end);
    Qt_neg_end = find(Qt_sign == 1, 1) + Qt_neg_start - 2;
    Qt_pos_start = Qt_neg_end + 1;
    Qt_sign = Qt_sign(Qt_pos_start - Qt_neg_start + 1: end);
    Qt_pos_end = find(Qt_sign ~= 1, 1) + Qt_pos_start - 2;
end
Qt_pos = [Qt_pos_start: Qt_pos_end]; % Indices for positive tricuspid flow
Qt_neg = [Qt_neg_start: Qt_neg_end]; % Indices for negative tricuspid flow
% Pulmonary
Qp_sign = sign(Q_p);
if(Qp_sign(1) <= 0)
    Qp_pos_start = find(Qp_sign == 1, 1);
    Qp_sign = Qp_sign(Qp_pos_start: end);
    Qp_pos_end = find(Qp_sign ~= 1, 1) + Qp_pos_start - 2;
    Qp_neg_start = Qp_pos_end + 1;
    Qp_sign = Qp_sign(Qp_neg_start - Qp_pos_start + 1: end);
    Qp_neg_end = find(Qp_sign == 1, 1) + Qp_neg_start - 2;
else
    Qp_neg_start = find(Qp_sign ~= 1, 1);
    Qp_sign = Qp_sign(Qp_neg_start: end);
    Qp_neg_end = find(Qp_sign == 1, 1) + Qp_neg_start - 2;
    Qp_pos_start = Qp_neg_end + 1;
    Qp_sign = Qp_sign(Qp_pos_start - Qp_neg_start + 1: end);
    Qp_pos_end = find(Qp_sign ~= 1, 1) + Qp_pos_start - 2;
end
Qp_pos = [Qp_pos_start: Qp_pos_end]; % Indices for positive pulmonary flow
Qp_neg = [Qp_neg_start: Qp_neg_end]; % Indices for negative pulmonary flow

% STENOSIS GRADING
% MS

MPG_m = trapz(t(Qm_pos), P_LA(Qm_pos) - P_LV(Qm_pos)) ./ (t(Qm_pos_end) - t(Qm_pos_start));

% AS
MPG_a = trapz(t(Qa_pos), P_LV(Qa_pos) - P_SA(Qa_pos)) / (t(Qa_pos_end) - t(Qa_pos_start));

% TS
MPG_t = trapz(t(Qt_pos), P_RA(Qt_pos) - P_RV(Qt_pos)) / (t(Qt_pos_end) - t(Qt_pos_start));

% PS
MPG_p = trapz(t(Qp_pos), P_RV(Qp_pos) - P_PA(Qp_pos)) / (t(Qp_pos_end) - t(Qp_pos_start));
Peak_PG_p = max(P_RV(Qp_pos) - P_PA(Qp_pos));


% REGURGITATION GRADING
% MR
RVol_m = -trapz(t(Qm_neg), Q_m(Qm_neg)); % Regurgitant volume over period where Qm is negative
SV_LA_pos = trapz(t(Qm_pos), Q_m(Qm_pos)); 
RF_m = 100 * (RVol_m / SV_LA_pos);

% AI
RVol_a = -trapz(t(Qa_neg), Q_a(Qa_neg));
SV_LV_pos = trapz(t(Qa_pos), Q_a(Qa_pos));
RF_a = 100 * (RVol_a / SV_LV_pos);

% TR
RVol_t = -trapz(t(Qt_neg), Q_t(Qt_neg));
SV_RA_pos = trapz(t(Qt_pos), Q_t(Qt_pos));
RF_t = 100 * (RVol_t / SV_RA_pos);

% PI
RVol_p = -trapz(t(Qp_neg), Q_p(Qp_neg));
SV_RV_pos = trapz(t(Qp_pos), Q_p(Qp_pos));
RF_p = 100 * (RVol_p / SV_RV_pos);




end_beat_i = find(t >= T, 1) - 1; % Index for end of one complete cardiac cycle

SV_LV_tot = max(V_LV) - min(V_LV); 
SV_RV_tot = max(V_RV) - min(V_RV);

EF_LV = SV_LV_tot / max(V_LV); 
EF_RV = SV_RV_tot / max(V_RV);

Q_m_maxima = findpeaks(Q_m(1:end_beat_i));

if(length(Q_m_maxima) == 2)
    E_A_ratio = Q_m_maxima(1) / Q_m_maxima(2);
elseif(length(Q_m_maxima) > 2) 
    while length(Q_m_maxima) > 2
        Q_m_maxima(Q_m_maxima == min(Q_m_maxima)) = [];
    end
    assert(length(Q_m_maxima) == 2, 'EAr bug 1 runSim');
    E_A_ratio = Q_m_maxima(1) / Q_m_maxima(2);
else
    E_A_ratio = -1; 
end


CO = SV_LV_pos * HR / 1000; % FIXME: Using CO with forward flow through LV 

MPAP = (1/3) * (max(P_PA)) + (2/3) * (min(P_PA));

[~, LVED_i] = max(V_LV); 
[~, RVED_i] = max(V_RV);
[~, LVES_i] = min(V_LV); 
[~, RVES_i] = min(V_RV);

LVIDd = -xm_LV(LVED_i) + xm_SEP(LVED_i) - 0.5 * (d_SW(LVED_i) + d_LW(LVED_i)); 
LVIDs = -xm_LV(LVES_i) + xm_SEP(LVES_i) - 0.5 * (d_SW(LVES_i) + d_LW(LVES_i));
RVIDd = xm_RV(RVED_i) - xm_SEP(RVED_i) - 0.5 * (d_SW(RVED_i) + d_RW(RVED_i)); 
RVIDs = xm_RV(RVES_i) - xm_SEP(RVES_i) - 0.5 * (d_SW(RVES_i) + d_RW(RVES_i));

%calculate systolic time and diastolic time (need time when valve opens)

Q_a_beat1=find(Q_a(1:end_beat_i));
sys_time=t(Q_a_beat1(end))-t(Q_a_beat1(1));

dias_time=T-sys_time;

sys_dias_time=sys_time/dias_time;



%% Calculating cost function metrics
o = struct(); % Relevant model outputs, to be compared to target values

o.SBP = max(P_SA);
o.DBP = min(P_SA);
o.LVEDV = max(V_LV);
o.LVESV = min(V_LV);
o.EF= (max(V_LV)-min(V_LV))/max(V_LV);
o.EAr = E_A_ratio;
o.LAVmax = max(V_LA);
o.LAVmin = min(V_LA);
o.RVEDV = max(V_RV);
o.RVESV = min(V_RV);
o.RAVmax = max(V_RA);
o.RAVmin = min(V_RA);
o.RAPmax = max(P_RA);
o.RAPmin = min(P_RA);
o.RAPmean= mean(P_RA);
o.PPAs = max(P_PA);
o.PPAd = min(P_PA);
o.PCWP = mean(P_PV);
o.CVP = mean(P_SV);
o.CO = CO;
o.Hed_LW = d_LW(LVED_i);
o.Hed_SW = d_SW(LVED_i);
o.Hed_RW = d_RW(LVED_i);
o.RVEDP = P_RV(LVED_i);

o.MPG_m = MPG_m;
o.MPG_a = MPG_a;
o.MPG_t = MPG_t;
o.Peak_PG_p = Peak_PG_p;
o.RF_m = RF_m;
o.RF_a = RF_a;
o.RF_t = RF_t;
o.RF_p = RF_p;
o.sys_dias_time=sys_dias_time;

%% grading the valvular stuff
% 1] must be continuous and directionally monotonous (no local minima)
% I propose setting category means and interpolating the grades
% One would need to explain this in the paper, but I think it is not as
% complex. See this example

% grade mean targets, note the extrapolated target of 30 for hypothetical
% grade 3
% gmt = [2, (10-5)/2+5, 12, 25]; 
% g = [0, 1, 2, 3]; % grades
% 
% figure(101);clf;hold on;
% % plot the grade means as cursors
% plot(repmat(gmt, [2, 1]), repmat([0;10], [1, length(gmt)]), '--', 'Linewidth', 1);
% set(gca,'ColorOrderIndex',1)
% % for each target grade plot the square costs in the range of model outputs
% for gd = 0:2 % gd as for grade-data
%     mpg = 0:0.1:25; % mean pressure gradient - model output
%     % linear extrapolate of the grades - grade:model
%     % this is our model output grade for further cost function
%     gm = interp1(gmt, g, mpg, 'linear', 'extrap'); 
%     cost = (gd-gm).^2;
%     plot(mpg, cost, 'Linewidth', 2);
% end
% xlabel('MPG')
% ylabel('Cost')
% title('Cost for each target grades')
% legend('Grade 0 target', 'Grade 1 target', 'Grade 2 target', 'hypothetical Grade 3 target'...
%     , 'Target 0 - costs', 'Target 1 - costs', 'Target 2 - costs')



g = [1, 2, 3, 4]; % grades: none, mild, moderate, severe. FIXME: will need to preprocess data to change the 0 grades to 1.
g_t = [1, 2]; % Tricuspid grades. There is no grading convention for tricuspid stenosis, so this is supposed to mean 1 for no stenosis and 2 for stenosis.
% FIXME: see if any of these choices are good or not or if it matters.
% Mitral
gmt_MR = [0, 20, 40, 60];
gmt_MS = [2.5, 3.75, 7.5, 12]; 
% Aortic
gmt_AI = [0, 20, 40, 60];
gmt_AS = [2.25, 10, 30, 50];
% Tricuspid
gmt_TR = [0, 25, 35, 45];
gmt_TS = [0.5, 6]; % 
% Pulmonary
gmt_PI = [0, 10, 30 ,50];
gmt_PS = [0.67, 22, 50, 78]; % Peak pressure gradient



o.MR = interp1(gmt_MR, g, RF_m, 'linear', 'extrap');
o.MS = interp1(gmt_MS, g, MPG_m, 'linear', 'extrap');
o.AI = interp1(gmt_AI, g, RF_a, 'linear', 'extrap');
o.AS = interp1(gmt_AS, g, MPG_a, 'linear', 'extrap');
o.TR = interp1(gmt_TR, g, RF_t, 'linear', 'extrap'); 
o.TS = interp1(gmt_TS, g_t, MPG_t, 'linear', 'extrap'); 
o.PI = interp1(gmt_PI, g, RF_p, 'linear', 'extrap');
o.PS = interp1(gmt_PS, g, MPG_p, 'linear', 'extrap');

%% weights
% weighted by maximal target in a category
w = struct();
wp = 1/120; % weights for pressures    
w.SBP = wp; w.DBP = wp;
w.RAPmax = wp; w.RAPmin = wp; w.RAPmean = wp;
w.PPAs = wp*10; w.PPAd = wp*10;
w.PCWP = wp*5; w.CVP = wp*10; w.RVEDP = wp*10;


wv = 1/160; % weights for volumes           
w.LVEDV = wv; w.LVESV = wv;
w.LAVmax = wv; w.LAVmin = wv;
w.RVEDV = wv; w.RVESV = wv;
w.RAVmax = wv; w.RAVmin = wv;

% flows
w.CO = 10/targets.CO;

% ratios
w.EAr = 1/targets_rest.EAr/5;

wt = 1; % weight of thicknesses
w.Hed_LW = wt;
w.Hed_SW = wt;
w.Hed_RW = wt;


w.sys_dias_time=10;
w.EF=50;

wg = 0; % Weight of grades (FIXME)
w.MR = wg;
w.MS = wg;
w.AI = wg;
w.AS = wg;
w.TR = wg;
w.TS = wg;
w.PI = wg;
w.PS = wg;
%%

% List fields in both structs
outputsfn = fieldnames(o); 
targetsfn = fieldnames(targets); 

% Assign missing fields to s
N = length(targetsfn);
cost = zeros(1, N);



if ~exist('printStats')
    printStats = false;
end
for i = 1:N

% normalization by category - e^2/maxtargets
     cost(i) = (o.(targetsfn{i}) - targets.(targetsfn{i}))^2 * w.(targetsfn{i});


    if printStats
        fprintf('%i) %s: %1.2f (%1.2f)\t %1.0f%% (%1.2f€)\n', i, targetsfn{i}, o.(targetsfn{i}), targets.(targetsfn{i}), o.(targetsfn{i})/targets.(targetsfn{i})*100 - 100, cost(i));
    end
end
if printStats
    fprintf('Total cost: %1.2f \n\n', sum(cost));
end
totalCost = sum(cost);