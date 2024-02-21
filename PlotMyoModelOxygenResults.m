%% Plotting Cummulative Results
%Run after completion of Oxygen Transport Model at all exercise levels
% Or use results from AllResults.mat


Plot_All_R=1;
Plot_All_O2=1;
if Plot_All_R==1
    points=1:45;
    points=points-0.5;
    Max_R= [MyoModel_Results_Rest.R11max, MyoModel_Results_MildE1.R11max, MyoModel_Results_MildE2.R11max, MyoModel_Results_MildE3.R11max, MyoModel_Results_MaxE.R11max, ... 
            MyoModel_Results_Rest.R12max, MyoModel_Results_MildE1.R12max, MyoModel_Results_MildE2.R12max, MyoModel_Results_MildE3.R12max, MyoModel_Results_MaxE.R12max, ...
            MyoModel_Results_Rest.R13max, MyoModel_Results_MildE1.R13max, MyoModel_Results_MildE2.R13max, MyoModel_Results_MildE3.R13max, MyoModel_Results_MaxE.R13max, ...
            MyoModel_Results_Rest.Rm1max, MyoModel_Results_MildE1.Rm1max, MyoModel_Results_MildE2.Rm1max, MyoModel_Results_MildE3.Rm1max, MyoModel_Results_MaxE.Rm1max, ... 
            MyoModel_Results_Rest.Rm2max, MyoModel_Results_MildE1.Rm2max, MyoModel_Results_MildE2.Rm2max, MyoModel_Results_MildE3.Rm2max, MyoModel_Results_MaxE.Rm2max, ...
            MyoModel_Results_Rest.Rm3max, MyoModel_Results_MildE1.Rm3max, MyoModel_Results_MildE2.Rm3max, MyoModel_Results_MildE3.Rm3max, MyoModel_Results_MaxE.Rm3max, ...
            MyoModel_Results_Rest.R21max, MyoModel_Results_MildE1.R21max, MyoModel_Results_MildE2.R21max, MyoModel_Results_MildE3.R21max, MyoModel_Results_MaxE.R21max, ... 
            MyoModel_Results_Rest.R22max, MyoModel_Results_MildE1.R22max, MyoModel_Results_MildE2.R22max, MyoModel_Results_MildE3.R22max, MyoModel_Results_MaxE.R22max, ...
            MyoModel_Results_Rest.R23max, MyoModel_Results_MildE1.R23max, MyoModel_Results_MildE2.R23max, MyoModel_Results_MildE3.R23max, MyoModel_Results_MaxE.R23max];
    
    Min_R= [MyoModel_Results_Rest.R11min, MyoModel_Results_MildE1.R11min, MyoModel_Results_MildE2.R11min, MyoModel_Results_MildE3.R11min, MyoModel_Results_MaxE.R11min, ...
            MyoModel_Results_Rest.R12min, MyoModel_Results_MildE1.R12min, MyoModel_Results_MildE2.R12min, MyoModel_Results_MildE3.R12min, MyoModel_Results_MaxE.R12min, ...
            MyoModel_Results_Rest.R13min, MyoModel_Results_MildE1.R13min, MyoModel_Results_MildE2.R13min, MyoModel_Results_MildE3.R13min, MyoModel_Results_MaxE.R13min, ...
            MyoModel_Results_Rest.Rm1min, MyoModel_Results_MildE1.Rm1min, MyoModel_Results_MildE2.Rm1min, MyoModel_Results_MildE3.Rm1min, MyoModel_Results_MaxE.Rm1min, ...
            MyoModel_Results_Rest.Rm2min, MyoModel_Results_MildE1.Rm2min, MyoModel_Results_MildE2.Rm2min, MyoModel_Results_MildE3.Rm2min, MyoModel_Results_MaxE.Rm2min, ...
            MyoModel_Results_Rest.Rm3min, MyoModel_Results_MildE1.Rm3min, MyoModel_Results_MildE2.Rm3min, MyoModel_Results_MildE3.Rm3min, MyoModel_Results_MaxE.Rm3min, ...
            MyoModel_Results_Rest.R21min, MyoModel_Results_MildE1.R21min, MyoModel_Results_MildE2.R21min, MyoModel_Results_MildE3.R21min, MyoModel_Results_MaxE.R21min, ...
            MyoModel_Results_Rest.R22min, MyoModel_Results_MildE1.R22min, MyoModel_Results_MildE2.R22min, MyoModel_Results_MildE3.R22min, MyoModel_Results_MaxE.R22min, ...
            MyoModel_Results_Rest.R23min, MyoModel_Results_MildE1.R23min, MyoModel_Results_MildE2.R23min, MyoModel_Results_MildE3.R23min, MyoModel_Results_MaxE.R23min];
    
    figure()
    plot([points; points], [Min_R; Max_R], 'LineWidth',10)
    colororder(["#648FFF";"#5A45B9";"#DC267F";"#FF7B2A";"#966904"])
    set(gca,"FontSize",15)
    title("Resistance Across the Coronary Perfusion LPM")
    ylabel("Resistance (mmHg*s/mL)")
end


if Plot_All_O2==1
    CA= 0.2/22.4;
    PO2=linspace(0,50,500001);
    tol=1e-5;
    K1=0.0368;
    K2=0.0002;
    K3=2.224;
    K4=0.1730;
    O2SAT=(K1.*PO2+2*K1.*K2.*PO2.^2+3*K1.*K2.*K3.*PO2.^3+4.*K1.*K2.*K3.*K4.*PO2.^4)./(4*(1+K1.*PO2+K1.*K2.*PO2.^2+K1.*K2.*K3.*PO2.^3+K1.*K2.*K3.*K4.*PO2.^4));

    O2SAT_EPI=[O2_Rest.c_21_avg, O2_MildE1.c_21_avg, O2_MildE2.c_21_avg, O2_MildE3.c_21_avg, O2_MaxE.c_21_avg]/CA;
    O2SAT_MID=[O2_Rest.c_22_avg, O2_MildE1.c_22_avg, O2_MildE2.c_22_avg, O2_MildE3.c_22_avg, O2_MaxE.c_22_avg]/CA;
    O2SAT_ENDO=[O2_Rest.c_23_avg, O2_MildE1.c_23_avg, O2_MildE2.c_23_avg, O2_MildE3.c_23_avg, O2_MaxE.c_23_avg]/CA;
    
    O2SAT_Grad_EPI=[O2_Rest_Grad.c_21_avg, O2_MildE1_Grad.c_21_avg, O2_MildE2_Grad.c_21_avg, O2_MildE3_Grad.c_21_avg, O2_MaxE_Grad.c_21_avg]/CA;
    O2SAT_Grad_MID=[O2_Rest_Grad.c_22_avg, O2_MildE1_Grad.c_22_avg, O2_MildE2_Grad.c_22_avg, O2_MildE3_Grad.c_22_avg, O2_MaxE_Grad.c_22_avg]/CA;
    O2SAT_Grad_ENDO=[O2_Rest_Grad.c_23_avg, O2_MildE1_Grad.c_23_avg, O2_MildE2_Grad.c_23_avg, O2_MildE3_Grad.c_23_avg, O2_MaxE_Grad.c_23_avg]/CA;
    
    O2SAT_Grad_Total=(O2SAT_Grad_EPI+O2SAT_Grad_MID+O2SAT_Grad_ENDO)/3;

    for i=1:5
        PO2_EPI(i)=mean(PO2(abs(O2SAT-O2SAT_EPI(i)) < tol));
        PO2_MID(i)=mean(PO2(abs(O2SAT-O2SAT_MID(i)) < tol));
        PO2_ENDO(i)=mean(PO2(abs(O2SAT-O2SAT_ENDO(i)) < tol));
        PO2_Grad_EPI(i)=mean(PO2(abs(O2SAT-O2SAT_Grad_EPI(i)) < tol));
        PO2_Grad_MID(i)=mean(PO2(abs(O2SAT-O2SAT_Grad_MID(i)) < tol));
        PO2_Grad_ENDO(i)=mean(PO2(abs(O2SAT-O2SAT_Grad_ENDO(i)) < tol));
    end


    O2_Ratio_EPI=[O2_Rest.C_ratio_epi, O2_MildE1.C_ratio_epi, O2_MildE2.C_ratio_epi, O2_MildE3.C_ratio_epi, O2_MaxE.C_ratio_epi]*100;
    O2_Ratio_MID=[O2_Rest.C_ratio_mid, O2_MildE1.C_ratio_mid, O2_MildE2.C_ratio_mid, O2_MildE3.C_ratio_mid, O2_MaxE.C_ratio_mid]*100;
    O2_Ratio_ENDO=[O2_Rest.C_ratio_endo, O2_MildE1.C_ratio_endo, O2_MildE2.C_ratio_endo, O2_MildE3.C_ratio_endo, O2_MaxE.C_ratio_endo]*100;
    
    O2_Ratio_Grad_EPI=[O2_Rest_Grad.C_ratio_epi, O2_MildE1_Grad.C_ratio_epi, O2_MildE2_Grad.C_ratio_epi, O2_MildE3_Grad.C_ratio_epi, O2_MaxE_Grad.C_ratio_epi]*100;
    O2_Ratio_Grad_MID=[O2_Rest_Grad.C_ratio_mid, O2_MildE1_Grad.C_ratio_mid, O2_MildE2_Grad.C_ratio_mid, O2_MildE3_Grad.C_ratio_mid, O2_MaxE_Grad.C_ratio_mid]*100;
    O2_Ratio_Grad_ENDO=[O2_Rest_Grad.C_ratio_endo, O2_MildE1_Grad.C_ratio_endo, O2_MildE2_Grad.C_ratio_endo, O2_MildE3_Grad.C_ratio_endo, O2_MaxE_Grad.C_ratio_endo]*100;
    
    HR=[64,90,120,150,180];

    figure()
    hold on
    plot(HR,O2SAT_ENDO, '-o','Color','r' ,'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName', 'ENDO')
    plot(HR,O2SAT_MID, '-o','Color', 'g' ,'MarkerFaceColor', 'g','LineWidth', 2, 'DisplayName', 'MID')
    plot(HR,O2SAT_EPI, '-o','Color', 'b' ,'MarkerFaceColor', 'b', 'LineWidth', 2, 'DisplayName', 'EPI')
    set(gca,'fontsize',15); box on; grid on;
    ylim([0,0.4])

    ylabel("Venous Total Oxygen Content (%)")
    title("Constant Transmural MVO2")


    figure()
    hold on
    plot(HR,O2SAT_Grad_ENDO, '-o','Color','r','MarkerFaceColor', 'r',  'LineWidth', 2, 'DisplayName', 'ENDO')
    plot(HR,O2SAT_Grad_MID, '-o','Color','g' ,'MarkerFaceColor', 'g', 'LineWidth', 2, 'DisplayName', 'MID')
    plot(HR,O2SAT_Grad_EPI, '-o','Color','b' ,'MarkerFaceColor', 'b', 'LineWidth', 2, 'DisplayName', 'EPI')
%     plot(HR,O2SAT_Grad_Total, '-*', 'LineWidth', 2, 'DisplayName', 'EPI')
    set(gca,'fontsize',15); box on; grid on;
    ylim([0,0.4])
    title("Graded Transmural MVO2")
    ylabel("Venous Total Oxygen Content (%)")

    figure()
    hold on
    plot(HR,PO2_ENDO, '-o','Color', 'r' ,'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName', 'ENDO')
    plot(HR,PO2_MID, '-o', 'Color','g' ,'MarkerFaceColor', 'g','LineWidth', 2, 'DisplayName', 'MID')
    plot(HR,PO2_EPI, '-o','Color', 'b' ,'MarkerFaceColor', 'b','LineWidth', 2, 'DisplayName', 'EPI')
    set(gca,'fontsize',15); box on; grid on;
    ylim([0,30])
    title("Constant Transmural MVO2")
    ylabel("Partial Pressure of Oxygen (mmHg)")

    figure()
    hold on
    plot(HR,PO2_Grad_ENDO, '-o','Color','r' ,'MarkerFaceColor', 'r', 'LineWidth', 2, 'DisplayName', 'ENDO')
    plot(HR,PO2_Grad_MID, '-o','Color','g' ,'MarkerFaceColor', 'g', 'LineWidth', 2, 'DisplayName', 'MID')
    plot(HR,PO2_Grad_EPI, '-o','Color','b' ,'MarkerFaceColor', 'b', 'LineWidth', 2, 'DisplayName', 'EPI')
    set(gca,'fontsize',15); box on; grid on;
    ylim([0,30])
    title("Graded MVO2")
    ylabel("Partial Pressure of Oxygen (mmHg)")

end