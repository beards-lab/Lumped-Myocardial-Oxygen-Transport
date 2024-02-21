% plot the figures
%Using Paul Tols Colorblind friendly palette + black

% #332288 (blue), EE6677 (red), 228833(green), CCBB44 (Yellow), 
% 66CCEE(cyan), AA337 (PURPLE), #BBBBBB (grey), #000000 (black)


%% Left Ventricle and Aortic Pressure
figure(); hold on
plot(t, P_LV, 'color', "#EE6677", "LineWidth", 3)
plot(t, P_SA, 'color', "k", "LineWidth", 3)
legend('P_{LV}', 'P_{Aorta}')
axis([0 2*T 0 200]); grid;
% set(gca,"FontSize",15)
set(gca,'FontSize',15,'LineWidth',1,'TickDir','both','TickLength',...
    [0.01 0.05]);

box on

%% Right ventricle and Pulmonary Arterial Pressure
figure(); hold on
plot(t, P_RV, 'color', "#66CCEE", "LineWidth", 3)
plot(t, P_PA, 'color', "k", "LineWidth", 3)
legend('P_{RV}', 'P_{Pulm. Arteries}')
axis([0 2*T 0 35]); grid;
% set(gca,"FontSize",15)
set(gca,'FontSize',15,'LineWidth',1,'TickDir','both','TickLength',...
    [0.01 0.05]);

box on

%% Right and Left PV Loop
figure(); hold on
plot(V_LV, P_LV, 'color', "#EE6677", "LineWidth", 3)
plot(V_RV, P_RV, 'color', "#66CCEE", "LineWidth", 3)
legend('P_{Left}', 'P_{Right}')
axis([0 200 0 200]); grid;
% set(gca,"FontSize",15)
set(gca,'FontSize',15,'LineWidth',1,'TickDir','both','TickLength',...
    [0.01 0.05]);

box on 


