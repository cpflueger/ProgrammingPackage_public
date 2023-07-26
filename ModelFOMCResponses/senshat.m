%% Sensitivity function from Campbell and Cochrane (1999)
%% INPUTS: 
% shat: log deviation from steady state surplus consumption ratio
% Sbar: steady-state surplus consumption ratio (in levels) 
%% Output:
% lambda(shat,Sbar) as defined in (21) (appendix)

function L = senshat(shat, Sbar)
    sbar = log(Sbar);
    Smax = sbar+0.5*(1-Sbar^2);
    L    = real((1/Sbar)*sqrt(1-2*shat)-1).*(shat<=Smax-sbar);
end

