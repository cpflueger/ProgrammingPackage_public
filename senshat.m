%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute sensitivity function for log surplus consumption ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS: 
% shat: log deviation from steady state surplus consumption ratio
% Sbar: steady-state surplus consumption ratio (in levels) 
%% Output:
% lambda(shat)

function L = senshat(shat, Sbar)
    sbar = log(Sbar);
    Smax = sbar+0.5*(1-Sbar^2);
    L    = real((1/Sbar)*sqrt(1-2*shat)-1).*(shat<=Smax-sbar);
end

