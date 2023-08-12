%%%%%%%%%% Minimum working example for programming package %%%%%%%%%%%%%%%%
% Autors: John Y. Campbell, Carolin Pflueger, and Luis M. Viceira
% Title: Macroeconomic Drivers of Bond and Equity Risks (JPE, 2020)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Observation: If you have any question about this package, please 
% send it to us as a comment in the issues part of the GitHub repository 
% or send it to the following emails: pflueger.carolin@gmail.com - l.yepezsa95@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all
%% Define the classes macro1, num1 and asset
macro1=macro_dyn;
num1=num_set; 
asset=asset_p;

%% Parameter settings that can speed up or slow down code
% Define a dummy variable to run the risk neutral part (=1 for risk neutral)
asset.risk_neutral_run  =   1;

% Parameter settings for numerical solution method, e.g. grid density,
% number of simulations etc. 
num1 = num1.parameters;

% Fix random number generator if desired for exact replicability
rng('default');

%% Input parameters for Model

% Preference parameters
macro1.theta0 = 0.9658;         % Peristence surplus consumption
macro1.theta1 = -0.0500;        % Dependence Output Gap          
macro1.theta2 = 0.0200;         % Dependence Lagged Output Gap
macro1.phi    = 0.9300;         % Consumption-output gap                       
macro1.gamma  = 2;              % Utility curvature
macro1.g      = 0.004725;       % Consumption growth 
macro1.rf     = 0.00235;        % Risk-free rate

% Leverage parameter
macro1.delta  = 0.5000; 

% Correlations of model shocks (v_t)
macro1.corr_vNew = [ 1.0000,    0.1100,    0.1100; 
                     0.1100,    1.0000,   -0.1100;
                     0.1100,   -0.1100,   1.0000];  

% Matrix needed to determine the solution for the dynamics of the model
macro1.P = [  0.8134,   -0.7781,   -0.1682; 
             -0.1111,    0.5556,   -0.5556;
             -0.1111,   -0.1111,    0.3333];

% Standard deviations of model shocks (v_t)
macro1.sigma_vNew = [0.0017    0.0017    0.0014]; 

% Compute implied Euler equation coefficients
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

%% Solve and simulate risk-neutral asset prices

% Solve and simulate risk neutral asset prices. Only executed if
% macro1.risk_neutral_run ==1
disp("Solve and simulate asset prices for period 1:")
asset = asset.risk_neutral_ap(macro1, num1); 

%% Solve and simulate full asset prices.

% Implements the value function iteration and simulations
disp('Computing prices')
asset = asset.computeFn21(num1,macro1);   

% Simulate path for macroeconomic dynamics and asset prices
disp('Simulate moments')
asset = asset.SimulateMoments(num1,macro1);

%% Table 2
disp("---------------------Table 2 - Stocks------------------------------")
disp("Period 1 Model Moments")
asset.stocks

%% Replicate Table 2

% STOCKS
% Equity Premium: Quarterly log return on average stock return in excess of the log 3-month Treasury bill
% plus one-half times the log excess return variance
equityPremium                 = 4*(mean(asset.simulated.reteq) + .5*std(asset.simulated.reteq).^2/100);

% Stocks Volatility: Log excess stock return standard deviation (in annualized percent)
vol_stocks                    = std(asset.simulated.reteq)*2;

% Sharpe Ratio: Ratio between the equity premium and stocks volatility
sharpeRatio_stocks            = equityPremium/vol_stocks;

% Mean of the price-dividend ratio
meanPDlev                     = exp(mean(asset.simulated.pdlev));

% Volatility of the price-dividend ratio
stdDP                         = std(asset.simulated.pdlev); 

% AR(1) coeff. pd: Persistence of the price-dividend ratio
rhoDP                         = corrcoef(asset.simulated.pdlev(2:end,:), asset.simulated.pdlev(1:end-1,:));

% Coefficient and R^2 of 1-YR Excess Returns on pd: Predictability of annual stock returns from the lagged 
% price-dividend ratio (in annualized percent)
ret1yr                        = conv(asset.simulated.reteq,ones(1,4),'valid')/100;
[re1_coef,~,~,~,R2_re1]       = regress(ret1yr, [ones(size(asset.simulated.pdlev(1:end-4))), ...
                                asset.simulated.pdlev(1:end-4)]);
coeffRegRetOnPD1y             = re1_coef(2);
R2RegRetOnPD1y                = R2_re1(1);

% Replicating Table 2
disp("-----------------------Replicate Table 2 - Stocks---------------------------------")
var_names2                    = ["equityPremium"; "vol"; "sharpeRatio"; "meanPDlev"; "stdDP"; "rhoDP"; ...
                                 "coeffRegRetOnPD1y"; "R2RegRetOnPD1y"];   
var_vals2                     = [equityPremium vol_stocks sharpeRatio_stocks meanPDlev stdDP rhoDP(1,2) ...
                                 coeffRegRetOnPD1y R2RegRetOnPD1y]' ;
for i=1:size(var_names2)
    disp(sprintf('%s: %.4f ',var_names2(i),var_vals2(i)));
end

%% Table 3
disp("----------------------Table 3 - Bonds------------------------------")
disp("Period 1 Model Moments")
asset.nominalBonds

%% Replicate Table 3

% BONDS
% Term Premium: The average log 5-year nominal bond return in excess of the log nominal 
% 3-month T-bill plus one-half times the log excess return variance 
termPremium                   = 4*(mean(asset.simulated.retnom) + .5*std(asset.simulated.retnom).^2/100);

% Volatility Excess Returns: Log excess bond return standard deviation (in annualized percent)
vol_nbond                     = std(asset.simulated.retnom)*2;

% Sharpe Ratio: Ratio between the term premium and bonds volatility
sharpeRatio_nbonds            = termPremium /vol_nbond;

% Yield Spread: The log 5-year bond yield minus the log nominal 3-month Treasury bill
spreadNom                     = asset.simulated.y5nom'-asset.simulated.rfr_nom;
meanLogYieldSpread            = mean(spreadNom);

% Volatility Yield Spread: The standard deviation of bond yield spreads
volLogYieldSpread             = std(spreadNom);

% AR(1) coeff. ys: Persistence of the yield spread
persistenceLogYieldSpread     = corrcoef(spreadNom(2:end), spreadNom(1:end-1));

% Coefficient and R^2 of 1-YR Excess Returns on log ys: Predictability of annual stock returns from the lagged 
% yield spread (in annualized percent)
ret1yrNom                     = asset.simulated.retnom(4:end)+asset.simulated.retnom(3:end-1)...
                                +asset.simulated.retnom(2:end-2)+asset.simulated.retnom(1:end-3);
ret1yrNom                     = ret1yrNom/100;                  
[ys1_coef,~,~,~,R2_ys1]       = regress(100*ret1yrNom, [ones(size(spreadNom(1:end-4)')), spreadNom(1:end-4)']);
coeffRegRetOnYS1y             = ys1_coef(2);
R2RegRetOnYS1y                = R2_ys1(1);

% Replicating Table 3
disp("------------------------Replicate Table 3 - Bonds---------------------------------")
var_names3                    = ["termPremium"; "vol"; "sharpeRatio"; "meanLogYieldSpread"; "volLogYieldSpread";...
                                 "persistenceLogYieldSpread"; "coeffRegRetOnYS1y"; "R2RegRetOnYS1y"];   
var_vals3                     = [termPremium vol_nbond sharpeRatio_nbonds meanLogYieldSpread volLogYieldSpread ...
                                 persistenceLogYieldSpread(1,2) coeffRegRetOnYS1y R2RegRetOnYS1y]' ;
for i=1:size(var_names3)
    disp(sprintf('%s: %.4f ',var_names3(i),var_vals3(i)));
end

%% Table 4
disp("---------------------Table 4 - Bonds and Stocks--------------------")
disp("Period 1 Model Moments")
asset.crossAsset

%% Replicate Table 4

% BOND-STOCK COMOVEMENT
% Correlation Bond and Stock Returns: The correlation of quarterly log bond excess returns 
% with log stock excess returns
bondstock_corr_temp           = corrcoef(asset.simulated.retnom, asset.simulated.reteq);
tipsstock_corr_temp           = corrcoef(asset.simulated.retreal, asset.simulated.reteq);
correlations                  = [bondstock_corr_temp; tipsstock_corr_temp];
corrNomStock                  = correlations(1,2);

% Beta Bond Returns on Stock Returns: The slope coefficient from regressing quarterly 
% log bond excess returns onto log stock excess returns
beta_temp                     = regress(asset.simulated.retnom, [ones(size(asset.simulated.retnom,1),1), ...
                                asset.simulated.reteq(:,1)]);
betaNom                       = beta_temp(2);

% Replicating Table 4
disp("------------------------Replicate Table 4 - Bonds and Stocks----------------------")
var_names4                    = ["corrNomStock"; "betaNom"];   
var_vals4                     = [corrNomStock betaNom]' ;
for i=1:size(var_names4)
    disp(sprintf('%s: %.4f ',var_names4(i),var_vals4(i)));
end
    
