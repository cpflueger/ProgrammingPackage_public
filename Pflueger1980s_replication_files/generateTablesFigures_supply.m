%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file generates the main tables and figures from the paper "Back to the 1980s or Not? %%
% The Drivers of Inflation and Real Risks in Treasury Bonds"                            %% 
%% cpflueger@uchicago.edu 1/20/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change directory to the specified path
%cd('C:\Users\cpflueger\Documents\GitHub\InflationExpectationsFormation\InflationExpectationsFormation');
cd('C:\Users\luisyepezsa\OneDrive - The University of Chicago\Documents\GitHub\ProgrammingPackage_public\Pflueger1980s_replication_files');
%cd('C:\Users\pflue\OneDrive\Documents\GitHub\InflationExpectationsFormation');

% Load empirical results data required for analysis
%load EmpiricalBase.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supply shock calibration (first-period) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear previous results
clear variables
close all

% Load supply shock calibration
load Calibrationgridsearch_supply2.mat

% Set initial macroeconomic moments based on the empirical first-period results
macro.moments=macro_base_emp.moments;

% Set adaptive inflation expectations and leverage parameters chosen to match Campbell-Shiller regressions
macro.zeta=0.6;
macro.delta=0.5;

% Set out-of-equilibrium shock volatilities as the supply shock calibration
macro.sigma_vec_out=macro.sigma_vec;

% Set bond maturity in quarters expected (1-year bond)
num.h=4;

% Set seed for replication
rng(0);

% Set to 1 to active the richer price-wage dynamics (Appendix G results)
price_wage_shock=0; 

% Set standard deviation of the price-wage shock
if price_wage_shock==1
    macro.sigmap=0.55/400;
else
    macro.sigmap=0;
end

% Compute implied Euler equation coefficients and fill in other parameters
macro = macro.update_params();

% Solves for the macro dynamics of the model 
macro = macro.ModelPQ82(num);

% Compute the rotated state vector 
macro = macro.ScaledStateVector();

% Get macro impulse responses
macro = macro.MacroIRF(num);

% Save the macro class results for use in figures or tables
macro_base=macro;
%macro_base_1=macro; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supply shock calibration with zeta=0 - this is needed to find the optimal zeta \in {0,0.6}%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set key calibration parameters according to the 2000s calibration, 
% while maintaining the remaining parameters from the supply calibration
macro.zeta=0;           
% macro.gamma_x=1/4;      
% macro.gamma_pi=1.1;     
% macro.rho_i=0.8;        
%macro.sigma_vec=[34.8100    0.0306    0.0306    0.0010]*10^(-6);   

% Compute implied Euler equation coefficients and fill in other parameters
macro = macro.update_params();

% Solves for the macro dynamics of the model 
macro = macro.ModelPQ82(num);

% Get macro impulse responses
macro = macro.ScaledStateVector();

% Get macro impulse responses
macro = macro.MacroIRF(num);

% Save the macro class results for use in figures or tables
macro_base2=macro;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demand shock calibration (second-period) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete previous macro class with the calibration
clear macro

% Load demand shock calibration
load Calibrationgridsearch_demand2.mat

% Set initial macroeconomic moments based on the empirical second-period results
macro.moments=macro_demand_emp.moments;

% Set out-of-equilibrium shock volatilities as the demand shock calibration
macro.sigma_vec_out=macro.sigma_vec;

% Set bond maturity in quarters expected (1-year bond)
num.h=4;

% Set seed for replication
rng(0);

% Set standard deviation of the price-wage shock
if price_wage_shock==1
    macro.sigmap=0.55/400;
else
    macro.sigmap=0;
end

% Compute implied Euler equation coefficients and fill in other parameters
macro = macro.update_params();

% Solves for the macro dynamics of the model 
macro = macro.ModelPQ82(num);

% Compute the rotated state vector 
macro = macro.ScaledStateVector();

% Get macro impulse responses
macro = macro.MacroIRF(num);

% Save the macro class results for use in figures or tables
macro_demand=macro;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demand calibration with zeta=0.6 - this is needed to find the optimal zeta \in {0,0.6}%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set demand calibration with backward-looking inflation expectations
macro=macro_demand;
macro.zeta=0.6;

% Compute implied Euler equation coefficients and fill in other parameters 
macro = macro.update_params();

% Solves for the macro dynamics of the model 
macro = macro.ModelPQ82(num);

% Compute the rotated state vector 
macro = macro.ScaledStateVector();

% Get macro impulse responses
macro = macro.MacroIRF(num);

% Save the macro class results for use in figures or tables
macro_demand2=macro;

% Save the macro class results with the supply calibration for use in figures or tables
macro=macro_base;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5: Plot Macro Impulse Responses  %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = figure('Position', get(0, 'Screensize'));
        pause(0.00001);
        frame_h = get(handle(gcf),'JavaFrame');
        set(frame_h,'Maximized',1); 
macro_alt=macro_demand;

subplot(4,3,1);
        p = plot(macro.Irf1.x*(1/100)*sqrt(1/macro.sigma_vec(1)), 'linewidth',2);
        set(p,'Color', [0,0,0])
        hold on
        plot(macro_alt.Irf1.x*(1/100)*sqrt(1/macro_alt.sigma_vec(1)), '--r', 'linewidth',2);
        ylim([-1,1.6])
        xlim([0,20])
        plot(0*macro.Irf1.x,'-k');
        set(gca, 'FontSize', 24)
        ylabel('Output Gap (x)','fontweight','normal','fontsize',17)
        title('Demand Shock','fontsize',24);

subplot(4,3,4);
        p = plot(macro.Irf1.pi*(1/100)*sqrt(1/macro.sigma_vec(1)), 'linewidth',2);
        set(p,'Color', [0,0,0])
        hold on
        plot(macro_alt.Irf1.pi*(1/100)*sqrt(1/macro_alt.sigma_vec(1)), '--r', 'linewidth',2);
        ylim([-0.2,2])
        xlim([0,20])
        plot(0*macro.Irf1.x,'-k');
        set(gca, 'FontSize', 24) % to increase ticks font size on axes
        ylabel('Inflation (\pi)','fontweight','normal','fontsize',17)


subplot(4,3,7);
        p = plot(macro.Irf1.i*(1/100)*sqrt(1/macro.sigma_vec(1)), 'linewidth',2);
        set(p,'Color', [0,0,0])
        hold on
        plot(macro_alt.Irf1.i*(1/100)*sqrt(1/macro_alt.sigma_vec(1)), '--r', 'linewidth',2);
        plot(0*macro.Irf1.pi,'-k');
        set(gca, 'FontSize', 24) % to increase ticks font size on axes
        ylim([-0.2,1.6])
        xlim([0,20])
        ylabel('Policy Rate (i)','fontweight','normal','fontsize',17)

subplot(4,3,10);
        p = plot(macro.Irf1.r*(1/100)*sqrt(1/macro.sigma_vec(1)), 'linewidth',2);
        set(p,'Color', [0,0,0])
        hold on
        plot(macro_alt.Irf1.r*(1/100)*sqrt(1/macro_alt.sigma_vec(1)), '--r', 'linewidth',2);
        plot(0*macro.Irf1.pi,'-k');
        set(gca, 'FontSize', 24) % to increase ticks font size on axes
        ylim([-0.2,1.5])
        xlim([0,20])
        ylabel('Real Rate (r)','fontweight','normal','fontsize',17)

subplot(4,3,2);
        p = plot(macro.Irf2.x*(1/400)*sqrt(1/macro.sigma_vec(2)), 'linewidth',2);
        set(p,'Color', [0,0,0])
        hold on
        plot(macro_alt.Irf2.x*(1/400)*sqrt(1/macro_alt.sigma_vec(2)), '--r', 'linewidth',2);
        ylim([-1,1.6])
        xlim([0,20])
        plot(0*macro.Irf1.x,'-k');
        set(gca, 'FontSize', 24)
        title('Supply Shock','fontsize',24);

subplot(4,3,5);
        p = plot(macro.Irf2.pi*(1/400)*sqrt(1/macro.sigma_vec(2)), 'linewidth',2);
        set(p,'Color', [0,0,0])
        hold on
        plot(macro_alt.Irf2.pi*(1/400)*sqrt(1/macro_alt.sigma_vec(2)), '--r', 'linewidth',2);
        ylim([-0.2,2])
        xlim([0,20])
        plot(0*macro.Irf2.x,'-k');
        set(gca, 'FontSize', 24) % to increase ticks font size on axes


subplot(4,3,8);
        p = plot(macro.Irf2.i*(1/400)*sqrt(1/macro.sigma_vec(2)), 'linewidth',2);
        set(p,'Color', [0,0,0])
        hold on
        plot(macro_alt.Irf2.i*(1/400)*sqrt(1/macro_alt.sigma_vec(2)), '--r', 'linewidth',2);
        plot(0*macro.Irf1.pi,'-k');
        set(gca, 'FontSize', 24) % to increase ticks font size on axes
        ylim([-0.2,1.6])
        xlim([0,20])

subplot(4,3,11);
        p = plot(macro.Irf2.r*(1/400)*sqrt(1/macro.sigma_vec(2)), 'linewidth',2);
        set(p,'Color', [0,0,0])
        hold on
        plot(macro_alt.Irf2.r*(1/400)*sqrt(1/macro_alt.sigma_vec(2)), '--r', 'linewidth',2);
        plot(0*macro.Irf1.pi,'-k');
        set(gca, 'FontSize', 24) % to increase ticks font size on axes
        ylim([-1.5,1])
        xlim([0,20])

subplot(4,3,3);
        p = plot(macro.Irf3.x*(1/400)*sqrt(1/macro.sigma_vec(3)), 'linewidth',2);
        set(p,'Color', [0,0,0])
        hold on
        plot(macro_alt.Irf3.x*(1/400)*sqrt(1/macro_alt.sigma_vec(3)), '--r', 'linewidth',2);
        ylim([-1,1.6])
        xlim([0,20])
        plot(0*macro.Irf1.x,'-k');
        set(gca, 'FontSize', 24)
        title('MP Shock','fontsize',24);

subplot(4,3,6);
        p = plot(macro.Irf3.pi*(1/400)*sqrt(1/macro.sigma_vec(3)), 'linewidth',2);
        set(p,'Color', [0,0,0])
        hold on
        plot(macro_alt.Irf3.pi*(1/400)*sqrt(1/macro_alt.sigma_vec(3)), '--r', 'linewidth',2);
        ylim([-0.2,2])
        xlim([0,20])
        plot(0*macro.Irf2.x,'-k');
        legend('1980s Calibration','2000s Calibration','fontsize',14)
        set(gca, 'FontSize', 24) % to increase ticks font size on axes


subplot(4,3,9);
        p = plot(macro.Irf3.i*(1/400)*sqrt(1/macro.sigma_vec(3)), 'linewidth',2);
        set(p,'Color', [0,0,0])
        hold on
        plot(macro_alt.Irf3.i*(1/400)*sqrt(1/macro_alt.sigma_vec(3)), '--r', 'linewidth',2);
        plot(0*macro.Irf1.pi,'-k');
        set(gca, 'FontSize', 24) % to increase ticks font size on axes
        ylim([-0.2,1.6])
        xlim([0,20])

 subplot(4,3,12);
        p = plot(macro.Irf3.r*(1/400)*sqrt(1/macro.sigma_vec(3)), 'linewidth',2);
        set(p,'Color', [0,0,0])
        hold on
        plot(macro_alt.Irf3.r*(1/400)*sqrt(1/macro_alt.sigma_vec(3)), '--r', 'linewidth',2);
        plot(0*macro.Irf1.pi,'-k');
        set(gca, 'FontSize', 24) % to increase ticks font size on axes
        ylim([-0.2,1.5])
        xlim([0,20])

% Save Figure 5
saveas(h, './figures/IRFmacro_supply_demand_units_Pflueger.png','png'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2: Plot Local Projections for Inflation, Output Gap, and Fed Funds Rate %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supply Calibration - Output Gap vs Inflation
macro=macro_base;
h=figure;
p  = plot((1:21), macro.moments.b,'--','linewidth',2); 
set(p,'Color', [0,0,0])
hold on 
plot((1:21), macro.coeffx(1:21), 'b', 'linewidth',2);
plot(0*(1:21),'-k');
plot((1:21), macro.moments.u, '--k', (1:21), macro.moments.d, '--k');
ylabel('Output Gap (x)', 'fontsize',20)
xlabel('Quarters after Inflation (\pi)', 'fontsize',20)
xlim([0,12])
ylim([-1.5,1.5])
hold off
legend('Data','Model','Location','SouthEast')
saveas(h,'./figures/LeadLag_x_pre_Pflueger.png','png');   
%% Supply Calibration - Output Gap vs Policy Rate
h=figure;
p  = plot((1:21), macro.moments.bix,'--','linewidth',2); 
set(p,'Color', [0,0,0])
hold on 
plot((1:21), macro.coeffix(1:21), 'b', 'linewidth',2);
plot(0*(1:21),'-k');
plot((1:21), macro.moments.uix, '--k', (1:21), macro.moments.dix, '--k');
xlim([0,12])
ylabel('Output Gap (x)', 'fontsize',20)
xlabel('Quarters after Policy Rate (i)', 'fontsize',20)
hold off
legend('Data','Model','Location','SouthEast')
saveas(h,'./figures/LeadLag_i_x_pre_Pflueger.png','png'); 
%% Supply Calibration - Policy Rate vs Inflation
h=figure;
p  = plot((1:21), macro.moments.bi,'--','linewidth',2); 
set(p,'Color', [0,0,0])
hold on 
plot((1:21), macro.coeffi(1:21), 'b', 'linewidth',2);
plot(0*(1:21),'-k');
plot((1:21), macro.moments.ui, '--k', (1:21), macro.moments.di, '--k');
xlim([0,12])
ylabel('Policy Rate (i)', 'fontsize',20)
xlabel('Quarters after Inflation (\pi)', 'fontsize',20)
hold off
legend('Data','Model','Location','SouthEast')
saveas(h,'./figures/LeadLag_i_pre_Pflueger.png','png'); 
%% Demand Calibration - Output Gap vs Inflation
h=figure;
p  = plot((1:21), macro_demand.moments.b,'--','linewidth',2); 
set(p,'Color', [0,0,0])
hold on 
plot((1:21), macro_demand.coeffx(1:21), 'b', 'linewidth',2);
plot(0*(1:21),'-k');
plot((1:21), macro_demand.moments.u, '--k', (1:21), macro_demand.moments.d, '--k');
xlim([0,12])
ylim([-1.5,1.5])
ylabel('Output Gap (x)', 'fontsize',20)
xlabel('Quarters after Inflation (\pi)', 'fontsize',20)
hold off
legend('Data','Model','Location','SouthEast')
saveas(h,'./figures/LeadLag_x_post_Pflueger.png','png');   
 
%% Demand Calibration - Output Gap vs Policy Rate
h=figure;
p  = plot((1:21), macro_demand.moments.bix,'--','linewidth',2); 
set(p,'Color', [0,0,0])
hold on 
plot((1:21), macro_demand.coeffix(1:21), 'b', 'linewidth',2);
plot(0*(1:21),'-k');
plot((1:21), macro_demand.moments.uix, '--k', (1:21), macro_demand.moments.dix, '--k');
xlim([0,12])
ylabel('Output Gap (x)', 'fontsize',20)
xlabel('Quarters after Policy Rate (i)', 'fontsize',20)
hold off
legend('Data','Model','Location','SouthEast')
saveas(h,'./figures/LeadLag_i_x_post_Pflueger.png','png'); 
%% Demand Calibration - Policy Rate vs Inflation
h=figure;
p  = plot((1:21), macro_demand.moments.bi,'--','linewidth',2); 
set(p,'Color', [0,0,0])
hold on 
plot((1:21), macro_demand.coeffi(1:21), 'b', 'linewidth',2);
plot(0*(1:21),'-k');
plot((1:21), macro_demand.moments.ui, '--k', (1:21), macro_demand.moments.di, '--k');
xlim([0,12])
ylabel('Policy Rate (i)', 'fontsize',20)
xlabel('Quarters after Inflation (\pi)', 'fontsize',20)
hold off
legend('Data','Model','Location','SouthEast')
saveas(h,'./figures/LeadLag_i_post_Pflueger.png','png'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure A2: Plot Local Projections for Inflation, Output Gap, and Fed Funds Rate with zeta=0 for supply calibration %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supply Calibration - Output Gap vs Inflation
macro=macro_base2;
h=figure;
p  = plot((1:21), macro.moments.b,'--','linewidth',2); 
set(p,'Color', [0,0,0])
hold on 
plot((1:21), macro.coeffx(1:21), 'b', 'linewidth',2);
plot(0*(1:21),'-k');
plot((1:21), macro.moments.u, '--k', (1:21), macro.moments.d, '--k');
ylabel('Output Gap (x)', 'fontsize',20)
xlabel('Quarters after Inflation (\pi)', 'fontsize',20)
xlim([0,12])
ylim([-1.5,1.5])
hold off
legend('Data','Model','Location','SouthEast')
saveas(h,'./figures/LeadLag_x_pre_zeta_Pflueger.png','png');   
%% Supply Calibration - Output Gap vs Policy Rate
h=figure;
p  = plot((1:21), macro.moments.bix,'--','linewidth',2); 
set(p,'Color', [0,0,0])
hold on 
plot((1:21), macro.coeffix(1:21), 'b', 'linewidth',2);
plot(0*(1:21),'-k');
plot((1:21), macro.moments.uix, '--k', (1:21), macro.moments.dix, '--k');
xlim([0,12])
ylabel('Output Gap (x)', 'fontsize',20)
xlabel('Quarters after Policy Rate (i)', 'fontsize',20)
hold off
legend('Data','Model','Location','SouthEast')
saveas(h,'./figures/LeadLag_i_x_pre_zeta_Pflueger.png','png'); 
%% Supply Calibration - Policy Rate vs Inflation
h=figure;
p  = plot((1:21), macro.moments.bi,'--','linewidth',2); 
set(p,'Color', [0,0,0])
hold on 
plot((1:21), macro.coeffi(1:21), 'b', 'linewidth',2);
plot(0*(1:21),'-k');
plot((1:21), macro.moments.ui, '--k', (1:21), macro.moments.di, '--k');
xlim([0,12])
ylabel('Policy Rate (i)', 'fontsize',20)
xlabel('Quarters after Inflation (\pi)', 'fontsize',20)
hold off
legend('Data','Model','Location','SouthEast')
saveas(h,'./figures/LeadLag_i_pre_zeta_Pflueger.png','png'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tables 1 and 2: Parameter calibrations and asset pricing moment results %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve supply Calibration, i.e. period 1
% Set macro results from supply calibration
macro=macro_base;

% Compute the rotated state vector 
macro = macro.ScaledStateVector; 

% Generate the grid for surplus consumption ratio shat
num = num.chooseshatgrid(macro);

% Generate the grid for value function iteration 
num = num.creategrid(macro);

% Set risk_neutral_run=0 to speed up code (=1 for risk neutral)
asset.risk_neutral_run = 1;

% Solve and simulate asset prices
asset = asset.SolveAll(macro, num);
%J-value including distance of model and data Campbell-Shiller regressions
J2_base=macro.J+100*((asset.nominalBonds.coeffRegRetOnYS1y-2.551955)/1.382841)^2;

% Save the num and asset class results with the supply calibration for use in figures or tables
num_base=num;
asset_base=asset;

%% Solve supply calibration with zeta=0
% Set macro results from supply calibration with zeta=0
macro=macro_base2;

% Compute the rotated state vector 
macro = macro.ScaledStateVector; 

% Generate the grid for surplus consumption ratio shat
num = num.chooseshatgrid(macro);

% Generate the grid for value function iteration 
num = num.creategrid(macro);

% Set risk_neutral_run=0 to speed up code (=1 for risk neutral)
asset.risk_neutral_run = 1;

% Solve and simulate asset prices
asset = asset.SolveAll(macro, num);
J2_base2=macro.J+100*((asset.nominalBonds.coeffRegRetOnYS1y-2.551955)/1.382841)^2;

%J2_base2>J2_base => set zet=0.6 for period 2 calibration
J2_base
J2_base2
% Save the num and asset class results with the supply calibration with zeta=0 for use in figures or tables
num_base2=num;
asset_base2=asset;

%% Solve demand calibration, i.e. period 2
% Set macro results from demand calibration
macro=macro_demand;

% Generate the grid for surplus consumption ratio shat
num = num.chooseshatgrid(macro);

% Generate the grid for value function iteration 
num = num.creategrid(macro);

% Set risk_neutral_run=0 to speed up code (=1 for risk neutral)
asset.risk_neutral_run = 1;

% Solve and simulate asset prices
asset = asset.SolveAll(macro, num);
J2_demand=macro.J+100*((asset.nominalBonds.coeffRegRetOnYS1y-0.8599105)/1.177813)^2;

% Save the num and asset class results with the demand calibration for use in figures or tables
asset_demand=asset;
num_demand=num;

%% Solve demand calibration with zeta=0.6
% Set macro results from demand calibration with zeta=0.6
macro=macro_demand2;

% Generate the grid for surplus consumption ratio shat
num = num.chooseshatgrid(macro);

% Generate the grid for value function iteration 
num = num.creategrid(macro);

% Set risk_neutral_run=0 to speed up code (=1 for risk neutral)
asset.risk_neutral_run = 1;

% Solve and simulate asset prices
asset = asset.SolveAll(macro, num);
J2_demand2=macro.J+100*((asset.nominalBonds.coeffRegRetOnYS1y-0.8599105)/1.177813)^2;

%J2_demand<J2_demand2 => set zet=0 for period 1 calibration
J2_demand
J2_demand2
% Save the num and asset class results with the demand calibration with zeta=0.6 for use in figures or tables
asset_demand2=asset;
num_demand2=num;

%% Table 2 (Asset Pricing Moments) period 1 
asset=asset_base;
macro=macro_base;
Table3=[asset.stocks.equityPremium, asset.stocks.vol, asset.stocks.sharpeRatio, asset.stocks.rhoDP, asset.stocks.coeffRegRetOnPD1y, asset.stocks.R2RegRetOnPD1y]';
Table3=[Table3; 0;  [asset.nominalBonds.expReturn, asset.nominalBonds.vol, asset.macroDynamics.y10nom_Vol, asset.nominalBonds.betaNom, asset.realBonds.betaRealStock, asset.nominalBonds.coeffRegRetOnYS1y, asset.nominalBonds.R2RegRetOnYS1y]';0;asset.macroDynamics.consGrowthVol;asset.macroDynamics.iChangeVol; asset.macroDynamics.Epi10Changes_Vol 
 ];
          
Table3 = [Table3];
Table3 = array2table(round(Table3,2));
Table3.Properties.RowNames = {'Equity Premium'...
'Equity Vol'...
'Equity SR'...
'AR(1) pd'...
'1 YR Excess Returns on pd'...
'1 YR Excess Returns on pd (R^2)'...
'-'...
'Bond Premium'...
'Return Vol.'...
'Quarterly Std. 10-Year Nominal Yields'...
'Nominal Bond-Stock Beta'...
'Real Bond-Stock Beta'...
'1 YR Excess Returns on Slope'...
'1 YR Excess Returns on Slope (R^2)'...
'--'...
'Std. Annual Cons. Growth'...
'Std Annual Change Fed Funds Rate'...
'Quarterly Std. 10-Year Inflation Forecast'...
};
Table3.Properties.VariableNames = {'Model'};

Table3supply=Table3

%% Table 1 (Parameters) period 1
params_calibrated=[macro.g*400, macro.gamma, macro.rf*400, macro.theta0^4, macro.theta1, 0, macro.gamma_pi, macro.gamma_x*4, macro.rho_i, 0, macro.kappa*4, macro.rho_pi, 0, macro.phi, macro.delta, 0, sqrt(macro.sigma_vec(1)).*100, sqrt(macro.sigma_vec(2)).*400, sqrt(macro.sigma_vec(3)).*400]';
params_implied=[macro.betaq^4, macro.Sbar, exp(macro.smax), macro.rhoxm, macro.rhoxp, macro.psi/4,macro.theta2]';

Table1 = (round(params_calibrated,2));

Table1supply=Table1;

%% Table 2 (Asset Pricing Moments) period 2
asset=asset_demand;
macro=macro_demand;
Table3=[asset.stocks.equityPremium, asset.stocks.vol, asset.stocks.sharpeRatio, asset.stocks.rhoDP, asset.stocks.coeffRegRetOnPD1y, asset.stocks.R2RegRetOnPD1y]';
Table3=[Table3; 0;  [asset.nominalBonds.expReturn, asset.nominalBonds.vol, asset.macroDynamics.y10nom_Vol, asset.nominalBonds.betaNom, asset.realBonds.betaRealStock, asset.nominalBonds.coeffRegRetOnYS1y, asset.nominalBonds.R2RegRetOnYS1y]';0;asset.macroDynamics.consGrowthVol;asset.macroDynamics.iChangeVol; asset.macroDynamics.Epi10Changes_Vol 
 ];
          
Table3 = array2table(round(Table3,2));
Table3.Properties.RowNames = {'Equity Premium'...
'Equity Vol'...
'Equity SR'...
'AR(1) pd'...
'1 YR Excess Returns on pd'...
'1 YR Excess Returns on pd (R^2)'...
'-'...
'Bond Premium'...
'Return Vol.'...
'Quarterly Std. 10-Year Nominal Yields'...
'Nominal Bond-Stock Beta'...
'Real Bond-Stock Beta'...
'1 YR Excess Returns on Slope'...
'1 YR Excess Returns on Slope (R^2)'...
'--'...
'Std. Annual Cons. Growth'...
'Std Annual Change Fed Funds Rate'...
'Quarterly Std. 10-Year Inflation Forecast'...
};
Table3.Properties.VariableNames = {'Model'};

Table3demand=Table3

%% Table 1 period 2
params_calibrated=[macro.g*400, macro.gamma, macro.rf*400, macro.theta0^4, macro.theta1, 0, macro.gamma_pi, macro.gamma_x*4, macro.rho_i, 0, macro.kappa*4, macro.rho_pi, 0, macro.phi, macro.delta, 0, sqrt(macro.sigma_vec(1)).*100, sqrt(macro.sigma_vec(2)).*400, sqrt(macro.sigma_vec(3)).*400]';
params_implied=[macro.betaq^4, macro.Sbar, exp(macro.smax), macro.rhoxm, macro.rhoxp, macro.psi/4,macro.theta2]';

Table1 = (round(params_calibrated,2));

Table1demand=Table1;

%% Table 1 with both calibrations
Table1combined=table([Table1supply, Table1demand]);

Table1combined.Properties.RowNames = {...
        '$g$',...
        '$\gamma$',...
        '$\bar r$',...
        '$\theta_0$',...
        '$\theta_1$',...
        '-',...
        '$\gamma_\pi$',...
        '$\gamma_x$',...
        '$\rho_i$',...
        '--',...
        '$\kappa$',...
        '$\rho_\pi$',...
        '----',...     
        '$\phi$',...
        '$\delta$',...
        '---',...
        '$\sigma_{x}$',...
        '$\sigma_{\pi}$',...
        '$\sigma_{i}$'
        };

%parameter standard errors are computed and reported in gridsearch_supply2
%(period 1) and gridsearch_demand (period 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Appendix Table A3: Inflation Forecast Error Regressions by Subperiod %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Supply Calibration result
macro_base.coeff_CG

% Demand Calibration result
macro_demand.coeff_CG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Asset price impulse responses to MP, Supply and Demand shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Asset prices IFRs for the Supply Calibration - MP Shock (Third Shock)
% Set the number of simulations
num_base.Nsim = 2; 

% Set the length of IRFs
num_base.Tirf = 150;

% Set the asset results for supply calibration
asset=asset_base;

% Set to 1 to compute IRF without background noise
asset.testirf=1;

% Set MP Shock (third shock)
asset.initialShockVec = [0,0,1/400,0];

% Compute simulations to build impulse responses for asset prices
asset =  asset.simulateStructuralIRF(macro_base,num_base);

% Save third shock structural IRF
asset.Irf3=asset.Irftemp;

% Save the asset class results with the supply calibration for use in figures or tables
asset_base=asset;

%% Compute Asset prices IFRs for the Demand Calibration - MP Shock (Third Shock)
% Set the number of simulations
num_demand.Nsim = 2; 

% Set the length of IRFs
num_demand.Tirf = 150;

% Set to 1 to compute IRF without background noise
asset_demand.testirf=1;

% Set MP Shock (third shock)
asset_demand.initialShockVec = [0,0,1/400,0];

% Compute simulations to build impulse responses for asset prices
asset_demand =  asset_demand.simulateStructuralIRF(macro_demand,num_demand);

% Save third shock structural IRF
asset_demand.Irf3=asset_demand.Irftemp;

%% Compute Asset prices IFRs for the Supply Calibration - Supply Shock (Second Shock)
% Set the asset results for supply calibration
asset=asset_base;

% Set Supply Shock (second shock)
asset.initialShockVec = [0,1/400,0,0];

% Compute simulations to build impulse responses for asset prices
asset =  asset.simulateStructuralIRF(macro_base,num_base);

% Save second shock structural IRF
asset.Irf2=asset.Irftemp;

% Save the asset class results with the supply calibration for use in figures or tables
asset_base=asset;

%% Compute Asset prices IFRs for the Demand Calibration - Supply Shock (Second Shock)
% Set Supply Shock (second shock)
asset_demand.initialShockVec = [0,1/400,0,0];

% Compute simulations to build impulse responses for asset prices
asset_demand =  asset_demand.simulateStructuralIRF(macro_demand,num_demand);

% Save second shock structural IRF
asset_demand.Irf2=asset_demand.Irftemp;

%% Compute Asset prices IFRs for the Supply Calibration - Demand Shock (First Shock)
% Set the asset results for supply calibration
asset=asset_base;

% Set Demand Shock (first shock)
asset.initialShockVec = [1/100,0,0,0];

% Compute simulations to build impulse responses for asset prices
asset =  asset.simulateStructuralIRF(macro_base,num_base);

% Save first shock structural IRF
asset.Irf1=asset.Irftemp;

% Save the asset class results with the supply calibration for use in figures or tables
asset_base=asset;

%% Compute Asset prices IFRs for the Demand Calibration - Demand Shock (First Shock)
% Set Demand Shock (first shock)
asset_demand.initialShockVec = [1/100,0,0,0];

% Compute simulations to build impulse responses for asset prices
asset_demand =  asset_demand.simulateStructuralIRF(macro_demand,num_demand);

% Save first shock structural IRF
asset_demand.Irf1=asset_demand.Irftemp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Appendix Figure A4: Yield Spread Impulse Response %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_IRF_CampbellShiller_2(asset_base, asset_demand, strcat('./figures/IRF_CS_combined_Pflueger'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 6: Model Asset Price Impulse Responses %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_IRF_Bonds_smallRP_real(asset_base, asset_demand, strcat('./figures/IRF_BondsRP_Pflueger'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Covariance and its decomposition into shocks - Appendix Table A1 and A2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Table A1 and A2 for period 1
%The tables below report the top 2x2 block of each panel in the Table
%totals are simply the sums across rows and columns, respectively
asset_base.crossAsset.cov_nom
asset_base.crossAsset.cov_real

%Table A2, Panels A and B for period 2
asset=asset_demand;
asset_demand.crossAsset.cov_nom
asset_demand.crossAsset.cov_real

% Save intermediate results of supply and demand calibrations for further calculations
save calibration_results.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 7:  Nominal Bond Betas by Prevalent vs. Realized Shocks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load supply and demand calibrations
load calibration_results.mat

% 1980s calibration
macro_base.sigma_vec_out=[];
macro_base.P_out=[];
macro_base.Q_out=[];
[Table1_base, Table3_base, asset_base, macro_base, num_base]=generateTable3(macro_base, num);

% 2000s prevalent, 2000s realized
macro_test=macro_base;
macro_test.sigma_vec_out=[];
macro_test.P_out=[];
macro_test.Q_out=[];
macro_test.sigma_vec=macro_demand.sigma_vec;
[Table1_sigma, Table3_sigma, asset_sigma, macro_sigma, num_sigma]=generateTable3(macro_test, num);

% 1980s prevalent, 2000s realized
macro_test=macro_base;
macro_test.sigma_vec_out=macro_demand.sigma_vec;
macro_test.P_out=[];
macro_test.Q_out=[];
[Table1_sigma_out, Table3_sigma_out, asset_sigma_out, macro_sigma_out, num_sigma_out]=generateTable3(macro_test, num);

% Save results for Figure 7
% Make sure that asset.risk_neutral_run = 1 in generateTable3. Otherwise,
% risk-neutral returns are not calculated
Counterfactuals_prevalent=[Table3_base, Table3_sigma, Table3_sigma_out]; 
save counterfactuals_out_of_equilibrium.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4:  Counterfactual Bond-Stock Betas  - Panel A %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part of the code computes the Figure A7 - Panel A if price_wage_shock==1
% If this portion of the code does not run, copy-paste the next section of
% code in the command window
% Load intermediate results of supply and demand calibrations
load calibration_results.mat

% Supply Calibration
macro_base.sigma_vec_out=[];
macro_base.P_out=[];
macro_base.Q_out=[];
[Table1_base, Table3_base, asset_base, macro_base, num_base]=generateTable3(macro_base, num);

% Supply Calibration changing all the MP parameters 
macro_base.sigma_vec_out=[];
macro_base.P_out=[];
macro_base.Q_out=[];
macro_base_MP=macro_base;
macro_base_MP.gamma_x=0.5*(macro_demand.gamma_x+macro_base.gamma_x);
macro_base_MP.gamma_pi=0.5*(macro_demand.gamma_pi+macro_base.gamma_pi);
macro_base_MP.rho_i=0.5*(macro_demand.rho_i+macro_base.rho_i);
[Table1_base_MP, Table3_base_MP, asset_base_MP, macro_base_MP, num_base_MP]=generateTable3(macro_base_MP, num);

% Supply Calibration changing the MP inertia 
macro_base_MP1=macro_base;
macro_base_MP1.rho_i=0.5*(macro_demand.rho_i+macro_base.rho_i);
[Table1_base_MP1, Table3_base_MP1, asset_base_MP1, macro_base_MP1, num_base_MP1]=generateTable3(macro_base_MP1, num);

% Supply Calibration changing the MP output/Inflation weights
macro_base_MP2=macro_base;
macro_base_MP2.gamma_x=0.5*(macro_demand.gamma_x+macro_base.gamma_x);
macro_base_MP2.gamma_pi=0.5*(macro_demand.gamma_pi+macro_base.gamma_pi);
[Table1_base_MP2, Table3_base_MP2, asset_base_MP2, macro_base_MP2, num_base_MP2]=generateTable3(macro_base_MP2, num);

% Supply Calibration changing shocks 
macro_base_sigma=macro_base;
macro_base_sigma.sigma_vec=(0.5*(sqrt(macro_base.sigma_vec)+sqrt(macro_demand.sigma_vec))).^2;
[Table1_base_sigma, Table3_base_sigma, asset_base_sigma, macro_base_sigma, num_base_sigma]=generateTable3(macro_base_sigma, num);

% Supply Calibration changing inflation expectations
macro_base_zeta=macro_base;
macro_base_zeta.zeta=0.5*(macro_demand.zeta+macro_base.zeta);
[Table1_base_zeta, Table3_base_zeta]=generateTable3(macro_base_zeta, num);

% Change the format of the previous table to fit the format to produce
% Figure 4 - Panel A - copy-paste into figures/Figure 4.xlsx
Table3_base_comparative=[Table3_base, Table3_base_sigma, Table3_base_MP, Table3_base_MP1, Table3_base_MP2, Table3_base_zeta];
Table3_base_comparative_new=[Table3_base_comparative(10:11,:)];

% Save the results for Figure 4 - Panel A
save comparative_statics_base.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4:  Counterfactual Bond-Stock Betas  - Panel B %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part of the code computes the Figure A7 - Panel B when sigmap=0.55

% Load the results for Figure 4 - Panel A
load comparative_statics_base.mat

% Demand Calibration
[Table1_demand, Table3_demand, asset_demand, macro_demand, num_demand]=generateTable3(macro_demand, num);

% Demand Calibration changing all the MP parameters 
macro_demand_MP=macro_demand;
macro_demand_MP.gamma_x=0.5*(macro_demand.gamma_x+macro_base.gamma_x);
macro_demand_MP.gamma_pi=0.5*(macro_demand.gamma_pi+macro_base.gamma_pi);
macro_demand_MP.rho_i=0.5*(macro_demand.rho_i+macro_base.rho_i);
[Table1_demand_MP, Table3_demand_MP, asset_demand_MP, macro_demand_MP, num_demand_MP]=generateTable3(macro_demand_MP, num);

% Demand Calibration changing the MP inertia 
macro_demand_MP1=macro_demand;
macro_demand_MP1.rho_i=0.5*(macro_demand.rho_i+macro_base.rho_i);
[Table1_demand_MP1, Table3_demand_MP1, asset_demand_MP1, macro_demand_MP1, num_demand_MP1]=generateTable3(macro_demand_MP1, num);

% Demand Calibration changing the MP output/Inflation weights
macro_demand_MP2=macro_demand;
macro_demand_MP2.gamma_x=0.5*(macro_demand.gamma_x+macro_base.gamma_x);
macro_demand_MP2.gamma_pi=0.5*(macro_demand.gamma_pi+macro_base.gamma_pi);
[Table1_demand_MP2, Table3_demand_MP2, asset_demand_MP2, macro_demand_MP2, num_demand_MP2]=generateTable3(macro_demand_MP2, num);

% Demand Calibration changing shocks 
macro_demand_sigma=macro_demand;
macro_demand_sigma.sigma_vec=(0.5*(sqrt(macro_base.sigma_vec)+sqrt(macro_demand.sigma_vec))).^2;
[Table1_demand_sigma, Table3_demand_sigma, asset_demand_sigma, macro_demand_sigma, num_demand_sigma]=generateTable3(macro_demand_sigma, num);

% Demand Calibration changing inflation expectations
macro_demand_zeta=macro_demand;
macro_demand_zeta.zeta=0.5*(macro_demand.zeta+macro_base.zeta);
[Table1_demand_zeta, Table3_demand_zeta, asset_demand_zeta, macro_demand_zeta, num_demand_zeta]=generateTable3(macro_demand_zeta, num);

% Change the format of the previous table to fit the format to produce Figure 4 - Panel B
Table3_demand_comparative=[Table3_demand, Table3_demand_sigma, Table3_demand_MP, Table3_demand_MP1, Table3_demand_MP2, Table3_demand_zeta];
Table3_demand_comparative_new=[Table3_demand_comparative(10:11,:)];

% Save the results for Figure 4 - Panel B
save comparative_statics_demand.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparative statics for Figure 3 (CS predictability) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load calibration_results.mat
zeta_vec=(0:0.2:0.8);
Table1_base_zeta=zeros(18,max(size(zeta_vec)));
Table1_demand_zeta=zeros(18,max(size(zeta_vec)));
Table3_base_zeta=zeros(22,max(size(zeta_vec)));
Table3_demand_zeta=zeros(22,max(size(zeta_vec)));

for j=1:max(size(zeta_vec))
    j
    macro=macro_base;
    macro.zeta=zeta_vec(j);
    [Table1_base_zeta(:,j), Table3_base_zeta(:,j)]=generateTable3(macro, num);
    
    macro=macro_demand;
    macro.zeta=zeta_vec(j);
    [Table1_demand_zeta(:,j), Table3_demand_zeta(:,j)]=generateTable3(macro, num);
    
end
%copy-paste into /figures/Figure3.xlsx
save CS_plot.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Robustness: Appendix Table A4 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load calibration_results.mat
[Table1_base, Table3_base]=generateTable3(macro_base, num);
[Table1_demand, Table3_demand]=generateTable3(macro_demand, num);

%set kappa=0.019
macro_base_r5=macro_base;
macro_base_r5.kappa=0.019/4;
[Table1_base_r5, Table3_base_r5]=generateTable3(macro_base_r5, num);

%set kappa=0.019
macro_demand_r5=macro_demand;
macro_demand_r5.kappa=0.019/4;
[Table1_demand_r5, Table3_demand_r5]=generateTable3(macro_demand_r5, num);

%set phi=1
macro_base_r1=macro_base;
macro_base_r1.phi=1;
[Table1_base_r1, Table3_base_r1]=generateTable3(macro_base_r1, num);

macro_demand_r1=macro_demand;
macro_demand_r1.phi=1;
[Table1_demand_r1, Table3_demand_r1]=generateTable3(macro_demand_r1, num);

%set zeta=0
macro_base_r2=macro_base;
macro_base_r2.zeta=0;
[Table1_base_r2, Table3_base_r2]=generateTable3(macro_base_r2, num);

macro_demand_r2=macro_demand;
macro_demand_r2.zeta=0.6;
[Table1_demand_r2, Table3_demand_r2]=generateTable3(macro_demand_r2, num);

%set gamma=1
macro_base_r3=macro_base;
macro_base_r3.gamma=1;
[Table1_base_r3, Table3_base_r3]=generateTable3(macro_base_r3, num);

macro_demand_r3=macro_demand;
macro_demand_r3.gamma=1;
[Table1_demand_r3, Table3_demand_r3]=generateTable3(macro_demand_r3, num);

%set gamma=1 and phi=1
macro_base_r4=macro_base;
macro_base_r4.gamma=1;
macro_base_r4.phi=1;
[Table1_base_r4, Table3_base_r4]=generateTable3(macro_base_r4, num);

macro_demand_r4=macro_demand;
macro_demand_r4.gamma=1;
macro_demand_r4.phi=1;
[Table1_demand_r4, Table3_demand_r4]=generateTable3(macro_demand_r4, num);

Table3_base_robust=[Table3_base, Table3_base_r1, Table3_base_r3, Table3_base_r4, Table3_base_r2, Table3_base_r5];
Table3_demand_robust=[Table3_demand, Table3_demand_r1, Table3_demand_r3, Table3_demand_r4, Table3_demand_r2, Table3_demand_r5];

save robustness.mat
