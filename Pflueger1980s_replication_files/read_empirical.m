%% load empirical moments
clear all
%load CalibrationBase.mat;
[data_pre,var]=xlsread('empirical_responses.xlsx','pre');
[data_post,var]=xlsread('empirical_responses.xlsx','post');
b_pre=data_pre(:,1);
u_pre=data_pre(:,2);
d_pre=data_pre(:,3);

bi_pre=data_pre(:,4);
ui_pre=data_pre(:,5);
di_pre=data_pre(:,6);

bix_pre=data_pre(:,7);
uix_pre=data_pre(:,8);
dix_pre=data_pre(:,9);

b_post=data_post(:,1);
u_post=data_post(:,2);
d_post=data_post(:,3);

bi_post=data_post(:,4);
ui_post=data_post(:,5);
di_post=data_post(:,6);

bix_post=data_post(:,7);
uix_post=data_post(:,8);
dix_post=data_post(:,9);

bwage_post=data_post(:,10);
uwage_post=data_post(:,11);
dwage_post=data_post(:,12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%empirical target moments with standard deviations
%%
macro_base.moments.std_c=1.154405;
macro_base.moments.std_c_se=1.154405/sqrt(2*87);
macro_base.moments.std_i=2.259651;
macro_base.moments.std_i_se=2.259651/sqrt(2*87);
macro_base.moments.std_infl10=.4681591;
macro_base.moments.std_infl10_se=.4681591/sqrt(2*57);
macro_base.moments.b=b_pre;
macro_base.moments.u=u_pre;
macro_base.moments.d=d_pre;
macro_base.moments.b_se=(macro_base.moments.d-macro_base.moments.b)/1.96;
macro_base.moments.bi=bi_pre;
macro_base.moments.ui=ui_pre;
macro_base.moments.di=di_pre;
macro_base.moments.bi_se=(macro_base.moments.di-macro_base.moments.bi)/1.96;
macro_base.moments.bix=bix_pre;
macro_base.moments.uix=uix_pre;
macro_base.moments.dix=dix_pre;
macro_base.moments.bix_se=(macro_base.moments.dix-macro_base.moments.bix)/1.96;
macro_base.moments.bwage=zeros(21,1);
macro_base.moments.uwage=zeros(21,1);
macro_base.moments.dwage=zeros(21,1);
macro_base.moments.bwage_se=ones(21,1);

macro_demand.moments.std_c=1.153939 ;
macro_demand.moments.std_c_se=macro_demand.moments.std_c/sqrt(2*75);
macro_demand.moments.std_i=1.40653;
macro_demand.moments.std_i_se=macro_demand.moments.std_i/sqrt(2*75);
macro_demand.moments.std_infl10=.0890392;
macro_demand.moments.std_infl10_se=.0890392 /sqrt(2*75);
macro_demand.moments.b=b_post;
macro_demand.moments.u=u_post;
macro_demand.moments.d=d_post;
macro_demand.moments.b_se=(macro_demand.moments.d-macro_demand.moments.b)/1.96;
macro_demand.moments.bi=bi_post;
macro_demand.moments.ui=ui_post;
macro_demand.moments.di=di_post;
macro_demand.moments.bi_se=(macro_demand.moments.di-macro_demand.moments.bi)/1.96;
macro_demand.moments.bix=bix_post;
macro_demand.moments.uix=uix_post;
macro_demand.moments.dix=dix_post;
macro_demand.moments.bix_se=(macro_demand.moments.dix-macro_demand.moments.bix)/1.96;
macro_demand.moments.bwage=bwage_post;
macro_demand.moments.uwage=uwage_post;
macro_demand.moments.dwage=dwage_post;
macro_demand.moments.bwage_se=(macro_demand.moments.dwage-macro_demand.moments.bwage)/1.96;
macro_empirical1=macro_base;
macro_empirical2=macro_demand;

[data_post,var]=xlsread('empirical_responses.xlsx','post_pandemic');
beta_nom_rolling=data_post(:,1);
beta_real_rolling=data_post(:,2);
date_rolling=var(2:end,1);

macro_base_emp=macro_base;
macro_demand_emp=macro_demand;
clear macro_base macro_demand
save EmpiricalBase.mat;