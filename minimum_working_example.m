%%%%%%%%%% Minimum working example for programming package %%%%%%%%%%%%%%%%
% Autors: John Y. Campbell, Carolin Pflueger, and Luis M. Viceira
% Title: Macroeconomic Drivers of Bond and Equity Risks (JPE, 2020)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables
close all
%% Define the class macro1, num1 and asset
macro1=macro_dyn;
num1=num_set; 
asset=asset_p;
%% Input Parameters for the class macro1, asset and num1
% Define a dummy variable to run the risk neutral part (=1 to risk neutral)
asset.risk_neutral_run = 1;
% macro1 class parameters
macro1.theta0 = 0.9658; 
macro1.theta1 = -0.0500;                  
macro1.theta2 = 0.0200; 
macro1.phi    = 0.9300;                                
macro1.gamma  = 2; 
macro1.g      = 0.004725;  
macro1.rf     = 0.00235; 
macro1.delta  = 0.5000; 
macro1.corr_vNew = [ 1.0000,    0.1100,    0.1100; 
                     0.1100,    1.0000,   -0.1100;
                     0.1100,   -0.1100,   1.0000];  
macro1.P = [  0.8134,   -0.7781,   -0.1682; 
             -0.1111,    0.5556,   -0.5556;
             -0.1111,   -0.1111,    0.3333];
macro1.sigma_vNew = [0.0017    0.0017    0.0014]; 
% Update the Euler equation parameters
macro1 = macro1.update_params; 
% Filling in some num1 parameters
num1 = num1.parameters;
%% Solve for macrodynamics
rng('default');
% Solves for the macro dynamics of the model 
macro1 = macro1.ModelPQ82(num1); 
% Compute the state vector 
macro1  = macro1.ScaledStateVector; 
%% Update some specifications of class num1
num1 = num1.update_num(macro1);
%% Solve and simulate risk-neutral asset prices
disp("Solve and simulate asset prices for period 1:")
rng('default'); 
% Solve and simulate risk neutral asset prices
asset = asset.risk_neutral_ap(macro1, num1);
%% Solve and simulate asset prices with risk premia
disp('Computing prices')
% Implements the value function iteration
asset = asset.computeFn21(num1,macro1);   
disp('Simulate moments')
asset = asset.SimulateMoments(num1,macro1);
%% Reproduce some series that help to replicate some moments
simulated = asset.simulated;
%% Table 2
disp("---------------------Table 2 - Stocks------------------------------")
disp("Period 1 Model Moments")
asset.stocks
%% Replicate Table 2
%equityPremium
equityPremium                 = 4*(mean(asset.simulated.reteq) + .5*std(asset.simulated.reteq).^2/100);
%vol
vol_stocks                    = std(asset.simulated.reteq)*2;
%sharpeRatio
sharpeRatio_stocks            = equityPremium/vol_stocks;
%meanPDlev
meanPDlev                     = exp(mean(asset.simulated.pdlev));
%stdDP
stdDP                         = std(asset.simulated.pdlev); 
%rhoDP
rhoDP                         = corrcoef(asset.simulated.pdlev(2:end,:), asset.simulated.pdlev(1:end-1,:));
%coeffRegRetOnPD1y and R2RegRetOnPD1y
ret1yr                        = conv(asset.simulated.reteq,ones(1,4),'valid')/100;
[re1_coef,~,~,~,R2_re1]       = regress(ret1yr, [ones(size(asset.simulated.pdlev(1:end-4))), ...
                                asset.simulated.pdlev(1:end-4)]);
coeffRegRetOnPD1y             = re1_coef(2);
R2RegRetOnPD1y                = R2_re1(1);
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
%termPremium
termPremium                   = 4*(mean(asset.simulated.retnom) + .5*std(asset.simulated.retnom).^2/100);
%vol
vol_nbond                     = std(asset.simulated.retnom)*2;
%sharpeRatio
sharpeRatio_nbonds            = termPremium /vol_nbond;
%meanLogYieldSpread
spreadNom                     = asset.simulated.y5nom'-asset.simulated.rfr_nom;
meanLogYieldSpread            = mean(spreadNom);
%volLogYieldSpread
volLogYieldSpread             = std(spreadNom);
%persistenceLogYieldSpread
persistenceLogYieldSpread     = corrcoef(spreadNom(2:end), spreadNom(1:end-1));
%coeffRegRetOnYS1y and R2RegRetOnYS1y
ret1yrNom                     = asset.simulated.retnom(4:end)+asset.simulated.retnom(3:end-1)...
                                +asset.simulated.retnom(2:end-2)+asset.simulated.retnom(1:end-3);
ret1yrNom                     = ret1yrNom/100;                  
[ys1_coef,~,~,~,R2_ys1]       = regress(100*ret1yrNom, [ones(size(spreadNom(1:end-4)')), spreadNom(1:end-4)']);
coeffRegRetOnYS1y             = ys1_coef(2);
R2RegRetOnYS1y                = R2_ys1(1);
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
%corrNomStock and betaNom 
bondstock_corr_temp           = corrcoef(asset.simulated.retnom, asset.simulated.reteq);
tipsstock_corr_temp           = corrcoef(asset.simulated.retreal, asset.simulated.reteq);
correlations                  = [bondstock_corr_temp; tipsstock_corr_temp];
corrNomStock                  = correlations(1,2);
beta_temp                     = regress(asset.simulated.retnom, [ones(size(asset.simulated.retnom,1),1), ...
                                asset.simulated.reteq(:,1)]);
betaNom                       = beta_temp(2);
disp("------------------------Replicate Table 4 - Bonds and Stocks----------------------")
var_names4                    = ["corrNomStock"; "betaNom"];   
var_vals4                     = [corrNomStock betaNom]' ;
for i=1:size(var_names4)
    disp(sprintf('%s: %.4f ',var_names4(i),var_vals4(i)));
end
    
