%function to solve out for main AP moments from macro, num, and asset
function [Table1, Table3, asset, macro, num]=generateTable3(macro, num)
asset=asset_p();
rng(0);
% Std of FOMC date monetary policy surprise in annualized percent
asset.initialShockVec  = [0 0 0.0652 0]; 

sigma_vec_out_reset=macro.sigma_vec_out;
P_out_reset=macro.P_out;
Q_out_reset=macro.Q_out;
if(abs(macro.sigma_vec_out)>0)
else
    macro.sigma_vec_out=macro.sigma_vec;
end

macro = macro.update_params();
macro = macro.ModelPQ82(num);
if(abs(macro.P_out)>0)
else
    macro.P_out=macro.P;
end
if(abs(macro.Q_out)>0)
else
    macro.Q_out=macro.Q;
end

macro = macro.ScaledStateVector();
num = num.chooseshatgrid(macro);
% Create main grid and numerical settings
num = num.creategrid(macro);

%set risk_neutral_run=0 to speed up code
%set risk_neutral_run=1 for risk premium decomposition 
asset.risk_neutral_run = 1;
asset = asset.SolveAll(macro, num);

%% generate Table 3
Table3=[asset.stocks.equityPremium, asset.stocks.vol, asset.stocks.sharpeRatio, asset.stocks.rhoDP, asset.stocks.coeffRegRetOnPD1y, asset.stocks.R2RegRetOnPD1y]';
Table3=[Table3; 0;  [asset.nominalBonds.expReturn, asset.nominalBonds.vol, asset.nominalBonds.betaNom, asset.realBonds.betaRealStock, asset.crossAsset.betaNom_rn,asset.crossAsset.betaReal_rn, asset.crossAsset.corrNomStock, asset.nominalBonds.coeffRegRetOnYS1y, asset.nominalBonds.R2RegRetOnYS1y]';0;asset.macroDynamics.consGrowthVol;asset.macroDynamics.iChangeVol; asset.macroDynamics.Epi10Changes_Vol; 0;0];
          
%% generate Table 1
params_calibrated=[macro.g*400, macro.gamma, macro.rf*400, macro.theta0^4, macro.theta1, macro.kappa*4, macro.phi, macro.delta, 0, macro.gamma_pi, macro.gamma_x*4, macro.rho_i, 0, sqrt(macro.sigma_vec(1)).*100, sqrt(macro.sigma_vec(2)).*400, sqrt(macro.sigma_vec(3)).*400, 0, macro.zeta ]';
Table1=round(params_calibrated,2);

%reset _out variables
macro.sigma_vec_out=sigma_vec_out_reset;
macro.P_out=P_out_reset;
macro.Q_out=Q_out_reset;
