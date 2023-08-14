%-------------------------------------------------------------------------%
% Minimum working example for programming package
% Autors: Carolin Pﬂueger and Gianluca Rinaldi
% Title: Why Does the Fed Move Markets so Much? 
% A Model of Monetary Policy and Time-Varying Risk Aversion (JFE, 2022)
%-------------------------------------------------------------------------%
% Main Observation: If you have any question about this package, please 
% send it to us as a comment in the issues part of the GitHub repository 
% or send it to the following emails: pflueger.carolin@gmail.com - l.yepezsa95@gmail.com
%-------------------------------------------------------------------------%
clear variables
close all
%% Define the classes macro1, num1 and asset
macro1  =   macro_dyn;
num1    =   num_set; 
asset   =   asset_p;

%% Parameter settings that can speed up or slow down code
% Define a dummy variable to run the risk neutral part (=1 for risk neutral)
asset.risk_neutral_run  =   1;

% Define a dummy variable to simulate the FOMC effects (=1 for FOMC effects)
asset.FOMC_run          =   1;

% Parameter settings for numerical solution method, e.g. grid density,
% number of simulations etc. 
num1 = num1.parameters;

% Fix random number generator if desired for exact replicability
rng('default');

%% Input parameters for New Keynesian model with FOMC surprises

% Monetary policy rule
macro1.gamma_x          =   0.5/4;          % MP Response to output (annualized percent)
macro1.gamma_pi         =   1.5;            % MP Response to inflation
macro1.rho_i            =   0.8;            % MP Persistence

% Preference parameters
macro1.theta0           =   0.9658;         % Peristence surplus consumption
macro1.theta1           =   -0.90;          % Backward looking habit
macro1.phi              =   0.9300;         % Consumption-output gap                                
macro1.gamma            =   2;              % Utility curvature
macro1.g                =   0.004725;       % Consumption growth      
macro1.rf               =   0.00235;        % Risk-free rate

% Phillips curve
macro1.kappa            =   0.0062/4;       % Slope of the Phillips curve
macro1.rho_pi           =   0.8;            % Backward-Looking Coefficient

% Leverage parameter
macro1.delta            =   0.6666; 

% Variance of quarterly MP shock in natural units. 
macro1.sigma_vec(3)     =   8.7767e-06;     

% Std of FOMC date monetary policy surprise in annualized percent
asset.initialShockVec  = [0 0 0.0652 0]; 

% Compute implied Euler equation coefficients and fill in other parameters
% such as the industry portfolio betas
macro1 = macro1.update_params;

%% Solve for macroeconomic dynamics of the form Y_t=PY_{t-1}+Qv_t 

% Solves for the macro dynamics of the model 
macro1 = macro1.ModelPQ82(num1);

%% Prepare for asset price value function iteration

% Compute the rotated state vector 
macro1  = macro1.ScaledStateVector; 

% Numerical settings for value function iterations, some of which depend on
% rotated grid
num1 = num1.update_num(macro1);

%% Solve and simulate risk-neutral asset prices. 

% Solve and simulate risk neutral asset prices. Only executed if
% macro1.risk_neutral_run ==1
disp("Solve and simulate asset prices")
asset = asset.risk_neutral_ap(macro1, num1); 

%% Solve and simulate full asset prices.

% Implements the value function iteration and simulations
disp('Computing prices')
asset = asset.computeFn21(num1,macro1);   

% Simulate path for macroeconomic dynamics and asset prices
disp('Simulate moments')
asset = asset.SimulateMoments(num1,macro1);

%% Structural Impulse responses. This computes the consumption moments in Table 3

% Macroeconomic impulse responses do not need simulation. Set num.Nsim=2
% and asset.testirf=1 for speed. 
num1.Nsim = 2; 
num1.Tirf = 150;
asset.testirf=1;

% Simulations to build impulse responses for asset prices
asset =  asset.simulateStructuralIRF(macro1,num1);

% Lag trough consumption response
lag_trough=find(asset.Irftemp.c==min(asset.Irftemp.c))-1;

%% Table 3: Quarterly moments

% Collect simulated model moments
Table3      =   [asset.stocks.equityPremium, asset.stocks.vol, asset.stocks.sharpeRatio, asset.stocks.rhoDP, asset.stocks.coeffRegRetOnPD1y, asset.stocks.R2RegRetOnPD1y]';
Table3      =   [Table3; 0;  [asset.nominalBonds.meanLogYieldSpread, asset.nominalBonds.vol, asset.nominalBonds.coeffRegRetOnYS1y, asset.nominalBonds.R2RegRetOnYS1y]';...
                0;asset.macroDynamics.consGrowthVol;asset.macroDynamics.iChangeVol];
Table3      =   [Table3; 0; min(asset.Irftemp.c)/max(asset.Irftemp.i); lag_trough];

% Data moments
Table3Data  = [7.89; 
              16.79; 
              0.47; 
              0.92;
              -0.38;
              0.23;
              0; 
              1.87;  
              9.35;
              2.69;
              0.14
              0;
              1.50;
              1.35;
              0;
              -0.7;
              4];

% Create Table 3 to compare data and simulated model moments
Table3 = [Table3, Table3Data];
Table3 = array2table(round(Table3,2));
Table3.Properties.RowNames = {'Equity Premium'...
'Equity Vol'...
'Equity SR'...
'AR(1) pd'...
'1 YR Excess Returns on pd'...
'1 YR Excess Returns on pd (R^2)'...
'-'...
'Yield Spread'...
'Return Vol.'...
'1-YR Excess REturns on Yield Spread'...
'1-YR Excess REturns on Yield Spread (R^2)'...
'--'...
'Std. Annual Cons. Growth'...
'Std Annual Change Fed Funds Rate'...
'---'...
'Trough Effect Consumption'...
'Lag Trough'
};
Table3.Properties.VariableNames = {'Model', 'Data'};
Table3

%% Table 4: Model Stock Returns on FOMC Dates

% Collect simulated model moments
Table4              =   [asset.AP_responses.bernankeKuttner(1), asset.AP_responses.bernankeKuttner_dummy(1,1), asset.AP_responses.bernankeKuttner(2),...
                         asset.AP_responses.bernankeKuttner_dummy(1,2), asset.AP_responses.bernankeKuttner(3), asset.AP_responses.bernankeKuttner_dummy(1,3);
                         0,	asset.AP_responses.bernankeKuttner_dummy(2,1), 0, asset.AP_responses.bernankeKuttner_dummy(2,2), 0, asset.AP_responses.bernankeKuttner_dummy(2,3)];

% Create Table 4 to show simulated model moments
Table4              =   array2table(round(Table4,2));
Table4.Properties.RowNames = {
 'FF_Shock'...
 'FF_Shock_x_(FF_Shock>0)'
 };
Table4.Properties.VariableNames = {'Overall'...
 'Overall_'...
 'RN'...
 'RN_'...
 'RP'...
 'RP_'
 };
Table4

%% This section shows how to replicate Table 3 from the simulated macro and asset pricing series directly

% STOCKS
% Equity Premium: Quarterly log return on average stock return in excess of the log 3-month Treasury bill
% plus one-half times the log excess return variance
equityPremium                 = 4*(mean(asset.simulated.reteq) + .5*std(asset.simulated.reteq).^2/100);

% Stocks Volatility: Log excess stock return standard deviation (in annualized percent)
vol_stocks                    = std(asset.simulated.reteq)*2;

% Sharpe Ratio: Ratio between the equity premium and stocks volatility
sharpeRatio_stocks            = equityPremium/vol_stocks;

% AR(1) coeff. pd: Persistence of the price-dividend ratio
rhoDP                         = corrcoef(asset.simulated.pdlev(2:end,:), asset.simulated.pdlev(1:end-1,:));

% Coefficient and R^2 of 1-YR Excess Returns on pd: Predictability of annual stock returns from the lagged 
% price-dividend ratio (in annualized percent)
ret1yr                        = conv(asset.simulated.reteq,ones(1,4),'valid')/100;
[re1_coef,~,~,~,R2_re1]       = regress(ret1yr, [ones(size(asset.simulated.pdlev(1:end-4))), ...
                                asset.simulated.pdlev(1:end-4)]);
coeffRegRetOnPD1y             = re1_coef(2);
R2RegRetOnPD1y                = R2_re1(1);


% 10-YEAR NOMINAL BONDS
% Yield Spread: The log 10-year bond yield minus the log nominal 3-month Treasury bill
spreadNom                     = asset.simulated.y10nom'-asset.simulated.rfr_nom;
meanLogYieldSpread            = mean(spreadNom);

% Volatility Excess Returns: Log excess bond return standard deviation (in annualized percent)
vol_nbond                     = std(asset.simulated.retnom)*2;

% Coefficient and R^2 of 1-YR Excess Returns on Yield Spread: Predictability of annual stock returns from 
% the lagged Yield Spread (in annualized percent)
ret1yrNom                     = asset.simulated.retnom(4:end)+asset.simulated.retnom(3:end-1)...
                                +asset.simulated.retnom(2:end-2)+asset.simulated.retnom(1:end-3);
ret1yrNom                     = ret1yrNom/100;                  
[ys1_coef,~,~,~,R2_ys1]       = regress(100*ret1yrNom, [ones(size(spreadNom(1:end-4)')), spreadNom(1:end-4)']);
coeffRegRetOnYS1y             = ys1_coef(2);
R2RegRetOnYS1y                = R2_ys1(1);


% MACROECONOMIC DYNAMICS
% Std. Annual Consumption Growth: Standard deviation of 4-quarter real consumption growth (in percent)
consGrowthVol                 = 100*std(asset.simulated.csim(5:end)-asset.simulated.csim(1:end-4));

% Std. Annual Change fed funds Rate: Standard deviation of 4-quarter change fed funds rate
iChangeVol                    = std(asset.simulated.rfr_nom(5:end)-asset.simulated.rfr_nom(1:end-4));


% Replica_Table3: Replicating Table 3
Replica_Table3  =   [equityPremium, vol_stocks, sharpeRatio_stocks, rhoDP(1,2), coeffRegRetOnPD1y, R2RegRetOnPD1y]';
Replica_Table3  =   [Replica_Table3; 0;  [meanLogYieldSpread, vol_nbond, coeffRegRetOnYS1y, R2RegRetOnYS1y]';0;consGrowthVol;iChangeVol];
Replica_Table3  =   [Replica_Table3; 0; min(asset.Irftemp.c)/max(asset.Irftemp.i); lag_trough];
Replica_Table3  =   [Replica_Table3, Table3Data];
Replica_Table3  =   array2table(round(Replica_Table3,2));

Replica_Table3.Properties.RowNames  =   Table3.Properties.RowNames;
Replica_Table3.Properties.VariableNames =   Table3.Properties.VariableNames;
Replica_Table3

%% This section shows how to replicate Table 4 from the simulated macro and asset pricing series directly

% FF Shock: Is the stock market response to monetary policy shocks in model-simulated data
bernankeKuttnerTemp1        =   regress(asset.simulated.reteq_FOMC,[ones(size(asset.simulated.rfr_nom_FOMC)); asset.simulated.rfr_nom_FOMC]');
bernankeKuttnerTemp2        =   regress(asset.simulated.reteq_rn_FOMC,[ones(size(asset.simulated.rfr_nom_FOMC)); asset.simulated.rfr_nom_FOMC]');
bernankeKuttner             =   [bernankeKuttnerTemp1(2), bernankeKuttnerTemp2(2), bernankeKuttnerTemp1(2)-bernankeKuttnerTemp2(2)];


% FF Shock × (FF Shock>0): Is the stock market response to monetary policy shocks in model-simulated data when the shock is positive
bernankeKuttnerTemp1_dummy  =   regress(asset.simulated.reteq_FOMC,[asset.simulated.rfr_nom_FOMC; asset.simulated.rfr_nom_FOMC.*(asset.simulated.rfr_nom_FOMC>0);...
                                (asset.simulated.rfr_nom_FOMC>0); ones(size(asset.simulated.rfr_nom_FOMC))]');
bernankeKuttnerTemp2_dummy  =   regress(asset.simulated.reteq_rn_FOMC,[ asset.simulated.rfr_nom_FOMC; asset.simulated.rfr_nom_FOMC.*(asset.simulated.rfr_nom_FOMC>0);...
                                (asset.simulated.rfr_nom_FOMC>0); ones(size(asset.simulated.rfr_nom_FOMC))]');
bernankeKuttner_dummy       =   [bernankeKuttnerTemp1_dummy, bernankeKuttnerTemp2_dummy, bernankeKuttnerTemp1_dummy-bernankeKuttnerTemp2_dummy];


% Replica_Table4: Replicating Table 4
Replica_Table4              =   [bernankeKuttner(1), bernankeKuttner_dummy(1,1), bernankeKuttner(2), bernankeKuttner_dummy(1,2), bernankeKuttner(3), bernankeKuttner_dummy(1,3);
                                 0,	bernankeKuttner_dummy(2,1), 0, bernankeKuttner_dummy(2,2), 0, bernankeKuttner_dummy(2,3)];

Replica_Table4              =   array2table(round(Replica_Table4,2));

Replica_Table4.Properties.RowNames      =   Table4.Properties.RowNames;
Replica_Table4.Properties.VariableNames =   Table4.Properties.VariableNames;
Replica_Table4

%% Figure 3: Asset price impulse responses with simulations
% Set the parameters for simulations as the number of simulations and 
% the length of IRFs
num1.Nsim = 5000; 
num1.Tirf = 150;
asset.testirf=0;

% Simulations to build impulse responses for asset prices
asset =  asset.simulateStructuralIRF(macro1,num1);

% Save third structural IRF
asset.Irf3=asset.Irftemp;

% Plotting function for Figure 3
plot_StructuralIRF_STMP(asset);
