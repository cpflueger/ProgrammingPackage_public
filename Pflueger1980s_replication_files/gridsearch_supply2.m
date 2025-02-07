%grid search for period 2 calibration
clear macro
%initial values
load Calibrationgridsearch_supply2_init.mat
%empirical target moments
load EmpiricalBase.mat
macro.moments=macro_base.moments;
rng(0)
macro_start=macro;

J=zeros(50,1);
random_variation=2*rand(max(size(J)),7);
random_variation(1,:)=ones(1,7);

macro_start.zeta=0;
sigma_vec=zeros(3,1);
sigma_vec(1)=macro_start.sigma_vec(1);
sigma_vec(2)=macro_start.sigma_vec(2);
sigma_vec(3)=macro_start.sigma_vec(3);

rho_i=macro_start.rho_i;
gamma_x=macro_start.gamma_x;
gamma_pi=macro_start.gamma_pi;
phi=macro_start.phi;

%loop over potential parameter values to find minimum difference between
%empirical and model macro moments
for i=1:max(size(J))
 macro.sigma_vec(1)=(round(sqrt(random_variation(i,1)*sigma_vec(1))*100,2)/100)^2;
 macro.sigma_vec(2)=(round(sqrt(random_variation(i,2)*sigma_vec(2))*400,2)/400)^2;
 macro.sigma_vec(3)=(round(sqrt(random_variation(i,3)*sigma_vec(3))*400,2)/400)^2;
 
 macro.gamma_x=min(round(gamma_x*random_variation(i,4)*4,1)/4,0.55);
 macro.gamma_pi=round(1.1+random_variation(i,5)*(gamma_pi-1.1),2);
 macro.rho_i=round(min(0.5+random_variation(i,6)*(rho_i-0.5),0.95),2);
 
macro = macro.update_params();
macro = macro.ModelPQ82(num);

J(i)=macro.J;
end
%% find smallest J-value. 
i=find(J==min(J));
%%
 macro.sigma_vec(1)=(round(sqrt(random_variation(i,1)*sigma_vec(1))*100,2)/100)^2;
 macro.sigma_vec(2)=(round(sqrt(random_variation(i,2)*sigma_vec(2))*400,2)/400)^2;
 macro.sigma_vec(3)=(round(sqrt(random_variation(i,3)*sigma_vec(3))*400,2)/400)^2;
 
 macro.gamma_x=min(round(gamma_x*random_variation(i,4)*4,1)/4,0.5);
 macro.gamma_pi=round(1.1+random_variation(i,5)*(gamma_pi-1.1),2);
 macro.rho_i=round(min(0.5+random_variation(i,6)*(rho_i-0.5),0.95),2);
 
 %% report optimal parameter values
sqrt([macro.sigma_vec(1:3)]).*[100,400,400]
[macro.gamma_pi, macro.gamma_x*4, macro.rho_i]
%% solve macro dynamics 
macro = macro.update_params();
macro = macro.ModelPQ82(num);
macro = macro.ScaledStateVector();
%% save
clear macro_base macro_empirical1 macro_empirical2 h
save Calibrationgridsearch_supply2.mat;
%% get standard errors
% simulate variance-covariance matrix of empirical moments from model - this works only if rng(0) is switched off in macro_dyn
load Calibrationgridsearch_supply2.mat;
differences_zscore_vec=zeros(13,100);
for i=1:1000
macro = macro.ModelPQ82(num);
differences_zscore_vec(:,i)=macro.differences_zscore;
end
Vhat=corrcoef(differences_zscore_vec');
save Calibrationgridsearch_supply2.mat;
%%
load Calibrationgridsearch_supply2.mat;
macro.P_out=[];
macro.Q_out=[];
macro.sigma_vec_out=[];

NH=20;
H_sum=zeros(6,13);
for j=1:NH
[SE_diag, H] = getParametersSE(macro, num);
H_sum=H_sum+H;
end
H=H_sum/NH;

%standard errors with non-diagonal variance-covariance matrix
Mhat=(H(:,1:12)*H(:,1:12)');
SE=sqrt(diag(inv(Mhat)*H(:,1:12)*Vhat(1:12,1:12)*H(:,1:12)'*inv(Mhat)'))';
SE(1)=4*SE(1);
SE_temp=SE;
SE(2)=SE_temp(1);
SE(1)=SE_temp(2);

[macro.gamma_pi, macro.gamma_x*4, macro.rho_i, sqrt(macro.sigma_vec(1))*100, sqrt(macro.sigma_vec(2))*400, sqrt(macro.sigma_vec(3))*400;
    round(SE,2)] 
save Calibrationgridsearch_supply2.mat;