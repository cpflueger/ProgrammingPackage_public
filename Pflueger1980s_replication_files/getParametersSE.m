%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimate standard errors for the macro parameters in Table 1 assuming diagonal variance-covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%INPUTS:   macro
%           num
%           Nsim = number of simulations used to estimate derivatives

%%OUTPUTS:  baseParams = estimated parameter values
%           Vhatp = variance-covariance matrix of estimated parameter
%           values
%           baseParams_scaled = estimated parameter values, scaled and
%           ordered as in Table 1
%           standardErrors_scaled = standard errors for baseParams_scaled

function [SE_temp, H] = getParametersSE(macro, num)
epsilon=0.05;
epsilon_x=0.005;
%% compute numerical Jacobian with respect to gamma_x

H=zeros(6,max(size(macro.differences_zscore)));

%scale epsilon for gamma_x by 0.25 to adjust for units
macro_epsilon=macro;
macro_epsilon.gamma_x=macro.gamma_x-0.25*epsilon;
macro_epsilon = macro_epsilon.update_params();
macro_epsilon = macro_epsilon.ModelPQ82(num);
differences_zscore_min=macro_epsilon.differences_zscore;

macro_epsilon.gamma_x=macro.gamma_x+0.25*epsilon;
macro_epsilon = macro_epsilon.update_params();
macro_epsilon = macro_epsilon.ModelPQ82(num);
differences_zscore_max=macro_epsilon.differences_zscore;

H(1,:)=0.5*(differences_zscore_max-differences_zscore_min)/(0.25*epsilon);
%% gamma_pi
macro_epsilon=macro;
macro_epsilon.gamma_pi=macro.gamma_pi-epsilon;
macro_epsilon = macro_epsilon.update_params();
macro_epsilon = macro_epsilon.ModelPQ82(num);
differences_zscore_min=macro_epsilon.differences_zscore;

macro_epsilon.gamma_pi=macro.gamma_pi+epsilon;
macro_epsilon = macro_epsilon.update_params();
macro_epsilon = macro_epsilon.ModelPQ82(num);
differences_zscore_max=macro_epsilon.differences_zscore;

H(2,:)=0.5*(differences_zscore_max-differences_zscore_min)/epsilon;
%% rho_i
macro_epsilon=macro;
macro_epsilon.rho_i=macro.rho_i-epsilon;
macro_epsilon = macro_epsilon.update_params();
macro_epsilon = macro_epsilon.ModelPQ82(num);
differences_zscore_min=macro_epsilon.differences_zscore;

macro_epsilon.rho_i=macro.rho_i+epsilon;
macro_epsilon = macro_epsilon.update_params();
macro_epsilon = macro_epsilon.ModelPQ82(num);
differences_zscore_max=macro_epsilon.differences_zscore;

H(3,:)=0.5*(differences_zscore_max-differences_zscore_min)/epsilon;

%% sigma_x
macro_epsilon=macro;
macro_epsilon.sigma_vec(1)=((sqrt(macro.sigma_vec(1))*100+epsilon_x)/100)^2;
macro_epsilon = macro_epsilon.update_params();
macro_epsilon = macro_epsilon.ModelPQ82(num);
differences_zscore_max=macro_epsilon.differences_zscore;

macro_epsilon.sigma_vec(1)=((sqrt(macro.sigma_vec(1))*100-epsilon_x)/100)^2;
macro_epsilon = macro_epsilon.update_params();
macro_epsilon = macro_epsilon.ModelPQ82(num);
differences_zscore_min=macro_epsilon.differences_zscore;

H(4,:)=(differences_zscore_max-differences_zscore_min)/(epsilon_x);

%% sigma_pi
macro_epsilon=macro;
macro_epsilon.sigma_vec(2)=((sqrt(macro.sigma_vec(2))*400+epsilon)/400)^2;
macro_epsilon = macro_epsilon.update_params();
macro_epsilon = macro_epsilon.ModelPQ82(num);
differences_zscore_max=macro_epsilon.differences_zscore;

macro_epsilon.sigma_vec(2)=((sqrt(macro.sigma_vec(2))*400-epsilon)/400)^2;
macro_epsilon = macro_epsilon.update_params();
macro_epsilon = macro_epsilon.ModelPQ82(num);
differences_zscore_min=macro_epsilon.differences_zscore;

H(5,:)=(differences_zscore_max-differences_zscore_min)/epsilon;
%% sigma_i
macro_epsilon=macro;
macro_epsilon.sigma_vec(3)=((sqrt(macro.sigma_vec(3))*400+epsilon)/400)^2;
macro_epsilon = macro_epsilon.update_params();
macro_epsilon = macro_epsilon.ModelPQ82(num);
differences_zscore_max=macro_epsilon.differences_zscore;

macro_epsilon.sigma_vec(3)=((sqrt(macro.sigma_vec(3))*400-epsilon)/400)^2;
macro_epsilon = macro_epsilon.update_params();
macro_epsilon = macro_epsilon.ModelPQ82(num);
differences_zscore_min=macro_epsilon.differences_zscore;

H(6,:)=(differences_zscore_max-differences_zscore_min)/epsilon;
%% standard errors
SE=sqrt(diag(inv(H*H')))';
SE(1)=4*SE(1);
SE_temp=SE;
SE(2)=SE_temp(1);
SE(1)=SE_temp(2);
SE_temp=SE;
end




