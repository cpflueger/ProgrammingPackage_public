%% This class contains all parameters and methods related to the macro dynamics of the model
classdef macro_dyn
    properties
      phi,                       % smoothing parameter for consumption 
      tau,                       % scaling parameter for consumption 
      g,                         % consumption growth rate
      delta,                     % leverage parameter relating the dividend growth to consumption growth 
      delta_betaport,            %XXXFL: leverage of beta-sorted portfolios
      rho_a,                     % dependence between expected productivity growth and lagged real rate
      theta0,                    % persistence of the surplus consumption ratio
      theta1,                    % dependence of the surplus consumption ratio on current output gap
      theta2,                    % dependence of the surplus consumption ratio on lagged output gap
      gamma,                     % parameter controlling utility curvature
      rf,                        % steady state real short-term interest rate at x_t = 0 
      alpha,                     % price-stickiness
      taucap,                    % capital share of production
      theta,                     % substitutability across varietal goods
      Frisch,                    % Frisch elasticity
      omega,                     % elasticity of marginal costs to own production
      sigma_vec,                 % variances of fundamental shocks, in particular, those are the calibrated
                                 % values for the variances of respectively, u^IS (set to 0), u^PC, u^MP and u^*
                                 % to get empirically relevant values in Table 4 of the paper, for example sigma^PC = 400*sqrt(sigma_vec(2))
      sigma_vec_out,             % out-of-equilibrium shock volatilities 
      gamma_x,                   % controls the reaction of the interest rate to the output gap
      gamma_pi,                  % governs the response of the interest rate to inflation 
                                 % relative to its target level
      gamma_a,                   % MP response to shock to potential output
      psi,                       % coefficient on (i_t - E(pi_t+1)) in IS equation
      rho_i,                     % governs the influence of past interest rates on current
                                 % interest rates (MP equation)
      rho_pi,                    % (1-rho_pi), rho_pi are, respectively, the coefficient on 
                                 % expected next period inflation and lagged inflation
                                 % in the PC equation
      zeta,                      % (1-zeta), zeta are, respectively, the weights on rational and adpative
                                 % components on inflation expectations
      rhoxm,                     % coefficient on lagged output gap in IS equation
      rhoxp,                     % coefficient on expected next period output gap in IS equation (tau/(tau*phi - theta1))
      kappa,                     % controls the sensitivity of inflation to the output gap
      P,                         % P, Q are the matrices determining the solution for the dynamics of 
                                 % Y_t, in the form Y_t = P*Y_t-1 + Q*u_t
      Q,                         %
      P_out,                     % out-of-equilibrium P 
      Q_out,                     % out-of- (asset pricing) equilibrium Q
      impact_negative,           % 1 if one of the diagonal entries of the solution Q is
                                 % negative, 0 otherwise
      number_stable,             % number of stable eigenvalues (modulus < 1)
      number_complex,            % number of complex eigenvalues
      eigenvalues,               % eigenvalues used to construct solution                
      eigenvalues_select,        % the three eigenvalues we actually select
%      explosive_prop,            % Proportion of simulations that result in explosive macro dynamics
      A,                         % Matrix rotating \hat Y_t so that shocks to Z tilde = A* \hat Y_t 
                                 % are independent standard normal and the first element is conditionally
                                 % perfectly correlated with consumption (Appendix A.2.3)
      Ainv,                      % Inverse of A  
      vast,                      % covariance of inflation target shocks with normalized shocks eps
      sigmaperp,                 % standard deviation of orthogonal component of u^* : refer to (61)-(65) in the appendix
      sigmap,                    % standard deviation of the new price-wage shock - LY 11/27/2024
      vecp,                      % covariance between sigmap_{t+1} and epsilon_{t+1} - LY 11/27/2024
      sigmapperp,                % standard deviation of orthogonal component of sigmap - LY 11/27/2024
      Ptilde,                    % Ptilde = A*P*A' 
      VarZtilde,                 % Unconditional variance of Z_t tilde = Ptilde* Z_{t-1} tilde + epsilon 
      StdZtilde,                 % Unconditional standard deviations of Z_t tilde (3X1 vector as we have 3 state variables) 
      sigmac,                    % Standard deviation of consumption shocks epsilon_c,t 
      betaq,                     % quarterly pure time discount rate 
      ImpliedParams,             % vector of parameters implied by the solution found
      Sbar,                      % steady state surplus-consumption ratio  
      sbar,                      % log of Sbar
      smax,                      % value of s (log surplus-consumption ratio) such that for larger s 
                                 % the sensitivity function equals zero (as in Campbell-Cochrane (1999))  
      Smax,                      % exp(smax)
      Sigmau,                    % variance covariance matrix of the fundamental shocks
      QM,                        % defined in equation (6) of the appendix
      empT,                      % sample size for this period
      data_5y_ffx_correlation,   % data correlation between 20-quarter average Fed Funds rate and output gap
      data_5y_infl_correlation,  % data correlation between 20-quarter average inflation rate and output gap
      data_ffx_correlation,      % data correlation between Fed Funds rate and output gap
      data_infl_correlation,     % data correlation between inflation and output gap
      J,                         % discrepancy between model and data IRFs similar to CEE (2005) but only for 1q,1y,3y,5y and 10y points
      solutionNumber,            % Identifies which solution has been selected by ModelPQ82
      Irf1,                      % IRF1 lists macro impulse responses to a one-standard deviation increase in shock 1
      Irf2,                      % IRF2 lists macro impulse responses to a one-standard deviation increase in shock 2
      Irf3,                      % IRF3 lists macro impulse responses to a one-standard deviation increase in shock 3
      Irf4,                      % IRF4 lists macro impulse responses to a one-standard deviation increase in shock 4
      coeffx,                    % regression coefficient a_{1,h} in x_{t+h}=a_{0,h}+a_{1,h} \pi_t+\varepsilon_{t+h}
      coeffxwage,                % regression coefficient a_{1,h} in x_{t+h}=a_{0,h}+a_{1,h} \piw_t+\varepsilon_{t+h}
      coeffi,                    % regression coefficient a_{1,h} in i_{t+h}=a_{0,h}+a_{1,h} \pi_t+\varepsilon_{t+h}
      coeffix,                   % regression coefficient a_{1,h} in x_{t+h}=a_{0,h}+a_{1,h} i_t+\varepsilon_{t+h}
      coeff_CG,                  % Coibion-Gorodnichenko coefficient
      rfrNomStd,                 %std annual changes in nominal risk-free rate
      consGrowthVol,             %std annual consumption growth
      pi1ChangeVol,              %std annual changes in one-year inflation expectations
      pi10ChangeVol,             %std annual changes in 10-year inflation expectations 
      moments,                   %empirical moments
      differences_zscore,        %z-score of differences between model and data moments
   end
   
   methods
       %% If one of the fundamental parameters is changed, updates the ones that depend on it
       function macro_dyn = update_params(macro_dyn)
           %coefficients in Euler equation add up to one with wage
           %inflation
           %macro_dyn.theta2=(macro_dyn.phi-macro_dyn.theta1)/(1+macro_dyn.psi*(1-macro_dyn.phi))-1;
           macro_dyn.theta2=(macro_dyn.phi-macro_dyn.theta1)-1;
           macro_dyn.psi    = 1/(macro_dyn.gamma*(macro_dyn.phi - macro_dyn.theta1));
           macro_dyn.rhoxp = 1/(macro_dyn.phi - macro_dyn.theta1);
           macro_dyn.rhoxm = macro_dyn.theta2/(macro_dyn.phi - macro_dyn.theta1);
       
           
           %back out discount rate and growth-adjusted discount rate from real
           %risk-free rate using relationship between equilbrium risk-free rate and
           %discount rate from Campbell and Cochrane (1999)
           logbetaq     = macro_dyn.gamma*macro_dyn.g-0.5*macro_dyn.gamma*(1-macro_dyn.theta0)-macro_dyn.rf;
           macro_dyn.betaq  = exp(logbetaq);
           
           macro_dyn.omega=(macro_dyn.taucap)/(1-macro_dyn.taucap)+1/macro_dyn.Frisch/(1-macro_dyn.taucap);

           %implied Phillips curve parameters
           betag=macro_dyn.betaq*exp(-(macro_dyn.gamma-1)*macro_dyn.g);
           %slope of PC
           %macro_dyn.kappa=((1-macro_dyn.alpha)/macro_dyn.alpha)*(1-betag*macro_dyn.alpha)/(1+betag)*(macro_dyn.omega/(1+macro_dyn.omega*macro_dyn.theta));
           macro_dyn.rho_pi=(1+macro_dyn.zeta*betag)/(1+betag); % Updated to inlcude hybrid inflation expectations: weighted average of rational and adaptive
           
           % LY 12/19/2024: Run this part only if sigmap is defined
           macro_dyn.sigma_vec(3)=macro_dyn.sigma_vec(3)+((1-macro_dyn.rho_i)^2)*(macro_dyn.gamma_pi^2)*(macro_dyn.sigmap^2);    % LY 11/27/2024: I updated the third component of the variance-covariance matrix of shocks according to Appendix equation (A148)
           if ~isempty(macro_dyn.sigma_vec_out)
            macro_dyn.sigma_vec_out(3)=macro_dyn.sigma_vec_out(3)+((1-macro_dyn.rho_i)^2)*(macro_dyn.gamma_pi^2)*(macro_dyn.sigmap^2);    % LY 11/27/2024: I made the same update for the out-of-equilibrium shock volatilities 
           end
       end
       
      %% Solves for the macro dynamics of the model
      function macro_dyn = ModelPQ82(macro_dyn, num_set)
        %% Set up matrices F,G,H,M, defining the equation to be solved    
        %wage price factor 
        factor_w=1/(1+macro_dyn.psi*(1-macro_dyn.phi));
        F= [macro_dyn.rhoxp*factor_w, macro_dyn.psi*factor_w       , 0;
            0              ,(1-macro_dyn.rho_pi) , 0;
            0              , 0                   , 0];

        G= [-1                                   , 0                                      ,-macro_dyn.psi*factor_w;
            macro_dyn.kappa                      , -1                                     ,0             ;
            (1-macro_dyn.rho_i)*macro_dyn.gamma_x, (1-macro_dyn.rho_i)*macro_dyn.gamma_pi ,-1           ];    
        

        H= [macro_dyn.rhoxm*factor_w, 0               , 0         ;
            0              , macro_dyn.rho_pi, 0               ;
            -(1-macro_dyn.phi)*(1-macro_dyn.rho_i)*macro_dyn.gamma_pi              , 0               , macro_dyn.rho_i];
        
        M= [factor_w, 0, 0, 0;
            0, 1, 0, 0;
            0, 0, 1, 1];
        
        %% Solve for P and Q
        m=3; % Number of state variables. In our case we have three: 
             % output gap, inflation and nominal interest rate gaps.

        % Selecting a solution is equivalent to picking three out of the 
        % six generalized eigenvalues, hence we impose a series of refinements as in
        % section A.2 of the appendix, recalling that eigenvalues with absolute
        % value <1 give stable solution

        % Construct matrices Xi and Delta as in Uhlig(1999) and compute the
        % generalized eigenvealues and eigenvectors 
        Xi=    [-G    , -H       ;
                eye(m), zeros(m)];

        Delta= [F       , zeros(m);
                zeros(m), eye(m) ];

        [X,e] = eig(Xi,Delta);
        evec  = diag(e);
               
        % Pick solution 
        [~,IX] = sort(abs(evec));  
        evec2  = evec(IX); 
        macro_dyn.eigenvalues = evec2;
        X2     = X(:,IX); 
        
        % Stable eigenvalues
        evec_stable = evec2(abs(evec2)<1);
        macro_dyn.number_stable  = size(evec_stable,1);
        %Count how many stable and complex eigenvalues there are
        evec_complex    = evec2(and(real(evec2)~=evec2,abs(evec2)<1));
        macro_dyn.number_complex = size(evec_complex,1); 

        % Implement the following with a switch
        % 3 stable evalues > unique solution
        % 4 stable evalues (2 complex) > 2 possible solutions
        % 4 stable evalues (all real) > 4 possible solutions
        % 5 stable evalues (2 complex) > 4 possible solutions
        % 5 stable evalues(all real) > 10 possible solutions

        switch macro_dyn.number_stable
            case 3
                select_mat = [true,true,true,false,false,false];
            case 4
                if macro_dyn.number_complex == 0
                    select_mat = [[true,true,true,false,false,false];
                              [true,true,false,true,false,false];
                              [true,false,true,true,false,false];
                              [false,true,true,true,false,false]];
                end
                if macro_dyn.number_complex == 2
                    select_complex = and(real(evec2)~=evec2,abs(evec2)<1);
                    select_complex = select_complex';
                    others = [1,1,1,1,0,0] - select_complex;
                    positions = find(others==1);
                    select1 = select_complex;
                    select1(positions(1)) = 1;
                    select2 = select_complex;
                    select2(positions(2)) = 1;
                    select_mat = [logical(select1);logical(select2)];
                end
                if macro_dyn.number_complex ==4
                    disp('No solutions!')
                end
            case 5
                if macro_dyn.number_complex == 0
                    select_mat = [[true,true,true,false,false,false];
                              [true,true,false,true,false,false];
                              [true,true,false,false,true,false];
                              [true,false,true,true,false,false];
                              [true,false,true,false,true,false];
                              [true,false,false,true,true,false];
                              [false,true,true,true,false,false];
                              [false,true,true,false,true,false];
                              [false,true,false,true,true,false];
                              [false,false,true,true,true,false]];
                end
                if macro_dyn.number_complex == 2
                    select_complex = and(real(evec2)~=evec2,abs(evec2)<1);
                    select_complex = select_complex';
                    others = [1,1,1,1,1,0] - select_complex;
                    positions = find(others==1);
                    select1 = select_complex;
                    select1(positions(1)) = 1;
                    select2 = select_complex;
                    select2(positions(2)) = 1;
                    select3 = select_complex;
                    select3(positions(3)) = 1;
                    select_mat = [logical(select1);logical(select2);logical(select3);logical(others)];
                    disp(select_mat)
                end
                if macro_dyn.number_complex ==4
                    select_complex = (and(real(evec2)~=evec2,abs(evec2)<1));
                    selected_complex = evec2(select_complex);
                    selected_complex = cplxpair(selected_complex);
                    other = [1,1,1,1,1,0] - select_complex';
                    select_complex1 = selected_complex(1:2);
                    select_complex2 = selected_complex(3:4);
                    select1 = other;
                    select1(evec2 == select_complex1(1)) = 1;
                    select1(evec2 == select_complex1(2)) = 1;  
                    select2 = other;
                    select2(evec2 == select_complex2(1)) = 1;
                    select2(evec2 == select_complex2(2)) = 1; 
                    select_mat = [logical(select1);logical(select2)];
                 end
        end
        
        Js         = zeros(size(select_mat,1),1);
        differences_zscore = zeros(size(select_mat,1),13);
        Ps         = zeros(size(select_mat,1),3,3);
        Qs         = zeros(size(select_mat,1),3,4);
        model_irfs = zeros(size(select_mat,1),num_set.maxT,9);
        explosive_p = zeros(size(select_mat,1),1);
        coeffx=zeros(size(select_mat,1),22);
        coeffxw=zeros(size(select_mat,1),22);
        coeffi=zeros(size(select_mat,1),22);
        coeffix=zeros(size(select_mat,1),22);
        CG_coeff=zeros(size(select_mat,1),1);
        for i = 1:size(select_mat,1)
            select = select_mat(i,:);
            % Issue warning if procedure does not pick exactly three eigenvalues 
            if(sum(select)~=3)
                warning('None or Multiple VAR(1) Solutions')
            end
            % Issue a warning if maximum selected eigenvalue has modulus greater than 1
            if max(abs(evec2(select)))>=1
                warning('Unstable Eigenvalue')
            end

            % Define matrices Omega and Lambda as in Uhlig(1999), using the  
            % eigenvectors corresponding to the selected eigenvalues
            Omega=X2(4:6,select); 
            Lambda=diag(evec2(select));

            macro_dyn.P = real(Omega*Lambda*inv(Omega)); %Compute solution for P 
            % Q has to satisfy Q = -[F*P +G]^(-1)*M so we can construct it as follows,
            % column by column
            macro_dyn.Q      = zeros(3,4);
            macro_dyn.Q(:,1) = -(F*macro_dyn.P+G)\M(:,1);
            macro_dyn.Q(:,2) = -(F*macro_dyn.P+G)\M(:,2);
            macro_dyn.Q(:,3) = -(F*macro_dyn.P+G)\M(:,3);
            macro_dyn.Q(:,4) = -(F*macro_dyn.P+G)\M(:,4);

            Ps(i,:,:)         = macro_dyn.P;
            Qs(i,:,:)         = macro_dyn.Q;
            
            % Calculate J value for this solution
             % Simulate irfs and calculate SMM objective function (J value) for this solution
                coeffx_sum = zeros(1,22);
                coeffxwage_sum = zeros(1,22);
                coeffi_sum = zeros(1,22);
                coeffix_sum = zeros(1,22);
                coeff_CG_sum = 0;
                rfrNomStd_sum= 0;
                consGrowthVol_sum=0;
                pi1ChangeVol_sum=0;
                pi10ChangeVol_sum=0;           
                
%                 explosive = 0;
                if(abs(macro_dyn.P_out)>0)
                    Pt = macro_dyn.P_out';
                else
                    %macro_dyn.P_out=macro_dyn.P;
                    Pt = macro_dyn.P'; % store the transpose of macro_dyn.P
                end
                
                if(abs(macro_dyn.Q_out)>0) 
                    Qt = macro_dyn.Q_out';
                else 
                    %macro_dyn.Q_out=macro_dyn.Q;
                    Qt = macro_dyn.Q'; % store the transpose of macro_dyn.Q
                end

                if(abs(macro_dyn.sigma_vec_out)>0)
                    sigma = diag((sqrt(macro_dyn.sigma_vec_out)).^2);
                    sqsigma = diag((sqrt(macro_dyn.sigma_vec_out)));
                else
                    %macro_dyn.sigma_vec_out=macro_dyn.sigma_vec;
                    sigma = diag((sqrt(macro_dyn.sigma_vec)).^2);
                    sqsigma = diag((sqrt(macro_dyn.sigma_vec)));
                end

                indexy = 1:3;
                
               
                %generate Nsim simulations from model and compute IRF from simulated data
                for j=1:num_set.Nsim
                    Y = zeros(macro_dyn.empT+2,4);
                    shk =  randn(macro_dyn.empT+1,size(sigma,1))*sqsigma;
                    shkQ = shk*Qt;
                    % Update state vector 
                    for t = 2:macro_dyn.empT+2   
                          Y(t,indexy) = Y(t-1,indexy)*Pt + shkQ(t-1,:); 
                    end
                    % Use cumsum to add the random walk component
                    Y(2:macro_dyn.empT+2,4)=cumsum(shk(:,4)); % We're now using this random shock component for productivity                
                    xt     = 100.*Y(:,1);
                    hatpit = 400.*Y(:,2);
                    hatit  = 400.*Y(:,3);
                    pistar = 100.*Y(:,4);   % Productivity (doesn't affect the other state variables)
                    pit    = hatpit;        % State variables already expressed in levels
                    it     = hatit;         % State variables already expressed in levels
                    %price inflation in annualized percent
                    pipt=pit;
                    pipt(2:end)=pipt(2:end)-4*(1-macro_dyn.phi)*xt(1:end-1);
                    %make sure price inflation gets plotted
                    piwt=pit;
                    pit=pipt;
                    %% lead-lag relationships inflation-output gap
                    coeffx_temp=zeros(1,22);
                    coeffxwage_temp=zeros(1,22);
                    coeffi_temp=zeros(1,22);
                    coeffix_temp=zeros(1,22);
                    for h=0:20
                    lhs    = xt(2+h:macro_dyn.empT);
                    rhs    = [ones(macro_dyn.empT-h-1,1), pit(2:macro_dyn.empT-h), pit(1:macro_dyn.empT-h-1)];
                    coeff_temp   = rhs\lhs;
                    coeffx_temp(h+1)=coeff_temp(2);
                    lhs    = xt(2+h:macro_dyn.empT);
                    rhs    = [ones(macro_dyn.empT-h-1,1), piwt(2:macro_dyn.empT-h), piwt(1:macro_dyn.empT-h-1)];
                    coeff_temp   = rhs\lhs;
                    coeffxwage_temp(h+1)=coeff_temp(2);
                    lhs    = it(2+h:macro_dyn.empT);
                    rhs    = [ones(macro_dyn.empT-h-1,1), pit(2:macro_dyn.empT-h), pit(1:macro_dyn.empT-h-1)];
                    coeff_temp   = rhs\lhs;
                    coeffi_temp(h+1)=coeff_temp(2);
                    lhs    = xt(2+h:macro_dyn.empT);
                    rhs    = [ones(macro_dyn.empT-h-1,1), it(2:macro_dyn.empT-h), it(1:macro_dyn.empT-h-1)];
                    coeff_temp   = rhs\lhs;
                    coeffix_temp(h+1)=coeff_temp(2);
                    end
                    %correlate long-term inflation and interest rates with
                    %output gap surprise
                    xshock=[0,[1,0,0]*shkQ'];
                    corrcoef_temp=corrcoef(xshock(1:end-39),conv(pit,ones(1,40),'valid'));
                    coeffx_temp(22)=corrcoef_temp(1,2);
                    corrcoef_temp=corrcoef(xshock(1:end-39),conv(it,ones(1,40),'valid'));
                    coeffi_temp(22)=corrcoef_temp(1,2);
                    corrcoef_temp=corrcoef(it(1:end-39),conv(pit,ones(1,40),'valid'));
                    coeffix_temp(22)=corrcoef_temp(1,2);
                     %subjective inflation expectations
                     pit_ma3= movmean(pit , 3 );
                     pit_ma4= movmean(pit , 4 );
                     Etilde3=macro_dyn.zeta*pit_ma3(1:end-2)+(1-macro_dyn.zeta)*400*(([0,1,0]*Pt'-(1-macro_dyn.phi)*[1,0,0])*(Y(3:end,1:3)*Pt^2)')';
                     Etilde4=macro_dyn.zeta*pit_ma4(2:end-1)+(1-macro_dyn.zeta)*400*(([0,1,0]*Pt'-(1-macro_dyn.phi)*[1,0,0])*(Y(3:end,1:3)*Pt^3)')';
                     revision=Etilde3(2:end)-Etilde4(1:end-1);
                     forecast_error=pit(6:end)-Etilde3(1:end-3);
                     lhs    = forecast_error(2:end);
                     rhs    = [ones(macro_dyn.empT-4,1), revision(1:end-3)];
                     coeff_temp   = rhs\lhs;
                     coeff_CG_temp=coeff_temp(2);
                     
                     %1-year inflation expectations
                     pi1 = macro_dyn.zeta*pit(2:end-1)'+(1-macro_dyn.zeta)*(1/4)*400*[0,1,0]*Pt'*(eye(3)-Pt'^4)*inv(eye(3)-Pt')*Y(3:end,1:3)';

                     pit_ma= movmean(pit, 40 );
                     %10-year inflation expectations
                     pi10 = macro_dyn.zeta*pit_ma(2:end-1)'+(1-macro_dyn.zeta)*(1/(4*10))*400*([0,1,0]*Pt'-(1-macro_dyn.phi)*[1,0,0])*(eye(3)-Pt'^(4*10))*inv(eye(3)-Pt')*Y(3:end,1:3)';
                     
                     %4-quarter consumption growth
                     Deltac_q = 100*[1,0,0]*Y(2:end,1:3)'-macro_dyn.phi*100*[1,0,0]*Y(1:end-1,1:3)';
                     Deltac_ann=Deltac_q(4:end)+Deltac_q(3:end-1)+Deltac_q(2:end-2)+Deltac_q(1:end-3);
                     
                     pi1Changes         = pi1(1:end-4) - pi1(5:end);
                     pi1ChangeVol_temp       = std(pi1Changes); 
                     
                     pi10Changes         = pi10(1:end-4) - pi10(5:end);
                     pi10ChangeVol_temp       = std(pi10Changes);    
               
                     % std Log consumption growth
                     consGrowthVol_temp     = std(Deltac_ann);
                     
                     rfr_nom              = 400*[0,0,1]*Y(:,1:3)';
                     rfrNomStd_temp       = std(rfr_nom(5:end)-rfr_nom(1:end-4));

                    coeffx_sum = coeffx_sum+coeffx_temp;
                    coeffxwage_sum = coeffxwage_sum+coeffxwage_temp;
                    coeffi_sum = coeffi_sum+coeffi_temp;
                    coeffix_sum = coeffix_sum+coeffix_temp;
                    coeff_CG_sum =coeff_CG_sum + coeff_CG_temp;
                    rfrNomStd_sum=rfrNomStd_sum+rfrNomStd_temp;
                    consGrowthVol_sum=consGrowthVol_sum+consGrowthVol_temp;
                    pi1ChangeVol_sum=pi1ChangeVol_sum+pi1ChangeVol_temp;
                    pi10ChangeVol_sum=pi10ChangeVol_sum+pi10ChangeVol_temp;             
                end   
                coeffx(i,:) = coeffx_sum/num_set.Nsim; 
                coeffxwage(i,:) = coeffxwage_sum/num_set.Nsim; 
                coeffi(i,:) = coeffi_sum/num_set.Nsim; 
                coeffix(i,:) = coeffix_sum/num_set.Nsim;
                coeff_CG(i) = coeff_CG_sum/num_set.Nsim;
                rfrNomStd(i) = rfrNomStd_sum/num_set.Nsim;
                consGrowthVol(i) = consGrowthVol_sum/num_set.Nsim;
                pi1ChangeVol(i) = pi1ChangeVol_sum/num_set.Nsim;
                pi10ChangeVol(i) = pi10ChangeVol_sum/num_set.Nsim;
                
               
                %objective function in terms of new moments    
                pricewage_model=coeffxwage(i,1)-coeffx(i,1);
                pricewage_data=macro_dyn.moments.bwage(1)-macro_dyn.moments.b(1);
                pricewage_se=sqrt(macro_dyn.moments.bwage_se(1)^2+macro_dyn.moments.b_se(1)^2);
                pricewage_zscore=(pricewage_model-pricewage_data)/pricewage_se;
                differences_new=[coeffx(i,2), coeffx(i,4), coeffx(i,10), coeffi(i,2), coeffi(i,4), coeffi(i,10), coeffix(i,2), coeffix(i,4), coeffix(i,10), consGrowthVol(i), rfrNomStd(i), pi10ChangeVol(i)];
                differences_new=differences_new-[macro_dyn.moments.b(2), macro_dyn.moments.b(4), macro_dyn.moments.b(10), macro_dyn.moments.bi(2), macro_dyn.moments.bi(4), macro_dyn.moments.bi(10), macro_dyn.moments.bix(2), macro_dyn.moments.bix(4), macro_dyn.moments.bix(10), macro_dyn.moments.std_c, macro_dyn.moments.std_i, macro_dyn.moments.std_infl10];
                se_new=[macro_dyn.moments.b_se(2), macro_dyn.moments.b_se(4), macro_dyn.moments.b_se(10), macro_dyn.moments.bi_se(2), macro_dyn.moments.bi_se(4), macro_dyn.moments.bi_se(10), macro_dyn.moments.bix_se(2), macro_dyn.moments.bix_se(4), macro_dyn.moments.bix_se(10), macro_dyn.moments.std_c_se, macro_dyn.moments.std_i_se, macro_dyn.moments.std_infl10_se];
                
                differences_new_zscore=differences_new./se_new;
                if abs(macro_dyn.moments.bwage(1))>10^(-10)
                    differences_new_zscore=[differences_new_zscore, pricewage_zscore];
                else
                    differences_new_zscore=[differences_new_zscore, 0];
                end
                Js(i) = differences_new_zscore*differences_new_zscore';
                differences_zscore(i,:)=differences_new_zscore;
                
                                %output gap shock loading onto u_t
                QMs(i,:) =[1,0,0]*macro_dyn.Q; 
%                 QM     = QMs(i,:);
%                 macro_dyn.Sigmau = diag(macro_dyn.sigma_vec_out);
        end
        %pick equilibrium with smallest distance between model and data impuilse responses 
        solution             = find(Js == min(Js));
        macro_dyn.eigenvalues_select = evec2(select_mat(solution,:));
        macro_dyn.P          = reshape(Ps(solution,:,:),3,3);
        macro_dyn.Q          = reshape(Qs(solution,:,:),3,4);
        macro_dyn.QM         = QMs(solution,:);
        macro_dyn.J          = Js(solution);
        macro_dyn.coeffx = coeffx(solution,:); 
        macro_dyn.coeffxwage = coeffxwage(solution,:);
        macro_dyn.coeffi = coeffi(solution,:); 
        macro_dyn.coeffix = coeffix(solution,:); 
        macro_dyn.coeff_CG = coeff_CG(solution);
        macro_dyn.rfrNomStd=rfrNomStd(solution);
        macro_dyn.consGrowthVol=consGrowthVol(solution);
        macro_dyn.pi1ChangeVol=pi1ChangeVol(solution);
        macro_dyn.pi10ChangeVol=pi10ChangeVol(solution);
        macro_dyn.differences_zscore=differences_zscore(solution,:);
        macro_dyn.solutionNumber = solution;
        %%
      end
     
      %% Compute the state vector Ztilde 
      function macro_dyn = ScaledStateVector(macro_dyn)
           
           % LY 12/19/2024: Run this part only if sigmap is not defined,
           % therefore we set sigmap=0 as default
           % if ~isfield(macro_dyn, 'sigmap') || isempty(macro_dyn.sigmap)
           %    macro_dyn.sigmap=0; 
           % end

           % Obtain the variance-covariance matrix of fundamental shocks
           macro_dyn.Sigmau = diag(macro_dyn.sigma_vec);

           % QM is loading of consumption shock onto u_t
           macro_dyn.QM     = macro_dyn.tau*[1,0,0]*macro_dyn.Q;
            
           % sigmac is the standard deviation of consumption shocks epsilon_c,t 
           macro_dyn.sigmac = sqrt(macro_dyn.QM*macro_dyn.Sigmau*macro_dyn.QM'); 

           % Implied steady-state surplus consumption ratio
           macro_dyn.Sbar   = macro_dyn.sigmac*sqrt(macro_dyn.gamma/(1-macro_dyn.theta0));
           %1+lambda(s) at s=steady-state
           lambda0          = 1/macro_dyn.Sbar;

           % Log surplus consumption steady state value
           macro_dyn.sbar   = log(macro_dyn.Sbar);

           % Maximum s_t
           macro_dyn.smax   = macro_dyn.sbar+0.5*(1-macro_dyn.Sbar^2);
           macro_dyn.Smax   = exp(macro_dyn.smax);

           % Back out betaq from risk-free rate and preferences parameters
           % as in Campbell and Cochrane (1999)
           logbetaq         = macro_dyn.gamma*macro_dyn.g-0.5*macro_dyn.gamma^2*macro_dyn.sigmac^2*lambda0^2-macro_dyn.rf;
           macro_dyn.betaq  = exp(logbetaq);
       
            % Save implied parameters
            macro_dyn.ImpliedParams = [macro_dyn.betaq^4, macro_dyn.rhoxm, macro_dyn.rhoxp, macro_dyn.psi, macro_dyn.Sbar,... 
                                       macro_dyn.smax, exp(macro_dyn.smax), macro_dyn.sigmac*200]';
            
            % variance-covariance matrix of Q u_t
            QSigma = macro_dyn.Q*macro_dyn.Sigmau*macro_dyn.Q';
            
            % consumption volatility
            macro_dyn.sigmac = sqrt(QSigma(1,1));
            
            % First row of A is proportional to first row of Q
            A1 = [1,0,0]/macro_dyn.sigmac;
            % pick A2 in the null space of A1*Q*Sigmau*Q'
            null1 = null(A1*macro_dyn.Q*macro_dyn.Sigmau*macro_dyn.Q');
            A2 = null1(:,1)';
            % scale so that A2*u has unit variance
            A2 = A2/sqrt(A2*macro_dyn.Q*macro_dyn.Sigmau*macro_dyn.Q'*A2');
            % pick A3 in null space of A1*Q*Sigmau*Q' and A2*Q*Sigmau*Q'
            null2 = null([A1*macro_dyn.Q*macro_dyn.Sigmau*macro_dyn.Q'; A2*macro_dyn.Q*macro_dyn.Sigmau*macro_dyn.Q']);
            A3 = null2';
            A3 = A3/sqrt(A3*macro_dyn.Q*macro_dyn.Sigmau*macro_dyn.Q'*A3');
            macro_dyn.A = [A1;A2;A3];
            
            % Compute A^{-1}
            macro_dyn.Ainv=inv(macro_dyn.A);

            % Compute loading of v_* onto normalized vector of shocks \eps_t. Corresponds
            % to vec_* in the appendix
            macro_dyn.vast=macro_dyn.A*macro_dyn.Q*[0,0,0,1]'*macro_dyn.sigma_vec(4);

            % sigmaperp = standard deviation of orthogonal component of v_* 
            macro_dyn.sigmaperp = sqrt(macro_dyn.sigma_vec(4)-macro_dyn.vast'*macro_dyn.vast);
            
            % LY 11/27/2024: Compute loading of v_p onto normalized vector of shocks \eps_t. Corresponds
            % to vec_p in the appendix - Equation (A150)
            macro_dyn.vecp=((1-macro_dyn.rho_i)*macro_dyn.gamma_pi*(macro_dyn.sigmap^2))*macro_dyn.A*macro_dyn.Q*[0,0,1,0]';
            
            % LY 11/27/2024: sigmapperp = standard deviation of orthogonal component of v_p - Equation (A151)
            macro_dyn.sigmapperp = sqrt((macro_dyn.sigmap^2)-macro_dyn.vecp'*macro_dyn.vecp);
            
            % dynamics for scaled state vector \tilde P
            macro_dyn.Ptilde                 = macro_dyn.A*macro_dyn.P/macro_dyn.A;
            
            % unconditional variance and standard deviation of \tilde Z
            % StdZtilde will be used determine the numerical grid for the value
            % function iteration
            macro_dyn.VarZtilde              = Uncon_var(macro_dyn.Ptilde, diag([1,1,1]));
            macro_dyn.StdZtilde              = sqrt(diag(macro_dyn.VarZtilde));
      end
      
            
      %compute macro impulse responses
      function macro_dyn = MacroIRF(macro_dyn,num_set)
          P = macro_dyn.P;
          Q = macro_dyn.Q;
          
            rng(0) %stop random generator so we can get exactly same results
            % Initialize matrices containing simulation results
            coeff = zeros(num_set.Nsim, 4,3);
            covar = zeros(num_set.Nsim,3,3);
            
            %choose which shock is the impulse
            for shocknum=1:4
          

                Y = zeros(macro_dyn.empT+1,4);
                for t = 2:macro_dyn.empT+1    % simulating data for chosen period/ empirical comparison

                    ut     = zeros(4,1);
                    if t==2
                        ut(shocknum)=sqrt(macro_dyn.sigma_vec(shocknum));
                    end
                    Y(t,1:3)   = transpose(P*Y(t-1,1:3)' + Q*ut); 
                    Y(t,4)     = Y(t-1,4) + ut(4); 
                end
                xt     = 100*Y(:,1);
                hatpit = 400*Y(:,2);
                hatit  = 400*Y(:,3);
                pistar = 100*Y(:,4); % Now productivity (doesn't affect the rest of state variables)
                pit = hatpit;
                it  = hatit;
                
%                 pipt=pit;
%                 pipt(2:end)=pipt(2:end)-4*(1-macro_dyn.phi)*xt(1:end-1);
%                 %make sure price inflation gets plotted
%                 piwt=pit;
%                 pit=pipt;
                
                
                e2 = [0,1,0];
                % 1 - quarter inflation expectations
                pi1qt = 400*e2*P*Y(:,1:3)'; 

                % 1, 2, 5, and 10-year inflation expectations
%                 pi1t = (1/4)*400*e2*P*(eye(3)-P^(4))*inv(eye(3)-P)*Y(:,1:3)'; 
%                 pi2 = (1/(2*4))*400*e2*P*(eye(3)-P^(2*4))*inv(eye(3)-P)*Y(:,1:3)';
%                 pi5 = (1/(5*4))*400*e2*P*(eye(3)-P^(5*4))*inv(eye(3)-P)*Y(:,1:3)';
%                 pi10= (1/(10*4))*400*e2*P*(eye(3)-P^(10*4))*inv(eye(3)-P)*Y(:,1:3)';

                rt = it-pi1qt'; % Real interest rate
                cgrowth = xt(2:end)-macro_dyn.phi*xt(1:end-1)+pistar(2:end)-pistar(1:end-1);
                ct = [0; cumsum(cgrowth/100)];
                Ct = exp(ct);

                %% save
                if shocknum==1
                macro_dyn.Irf1.x  = xt; 
                macro_dyn.Irf1.pi = pit;
                macro_dyn.Irf1.i  = it;
                macro_dyn.Irf1.r = rt;
                %macro_dyn.Irf1.pi1 = pi1t;
%                 macro_dyn.Irf1.pi2 = pi2;
%                 macro_dyn.Irf1.pi5 = pi5;
%                 macro_dyn.Irf1.pi10 = pi10;
                macro_dyn.Irf1.c = ct;
                macro_dyn.Irf1.C = Ct;
                end
                
                if shocknum==2
                macro_dyn.Irf2.x  = xt; 
                macro_dyn.Irf2.pi = pit;
                macro_dyn.Irf2.i  = it;
                macro_dyn.Irf2.r = rt;
                %macro_dyn.Irf2.pi1 = pi1t;
%                 macro_dyn.Irf2.pi2 = pi2;
%                 macro_dyn.Irf2.pi5 = pi5;
%                 macro_dyn.Irf2.pi10 = pi10;
                macro_dyn.Irf2.c = ct;
                macro_dyn.Irf2.C = Ct;
                end
                                
                if shocknum==3
                macro_dyn.Irf3.x  = xt; 
                macro_dyn.Irf3.pi = pit;
                macro_dyn.Irf3.i  = it;
                macro_dyn.Irf3.r = rt;
                %macro_dyn.Irf3.pi1 = pi1t;
%                 macro_dyn.Irf3.pi2 = pi2;
%                 macro_dyn.Irf3.pi5 = pi5;
%                 macro_dyn.Irf3.pi10 = pi10;
                macro_dyn.Irf3.c = ct;
                macro_dyn.Irf3.C = Ct;
                end
                
                if shocknum==4
                macro_dyn.Irf4.x  = xt; 
                macro_dyn.Irf4.pi = pit;
                macro_dyn.Irf4.i  = it;
                macro_dyn.Irf4.r = rt;
                %macro_dyn.Irf4.pi1 = pi1t;
%                 macro_dyn.Irf4.pi2 = pi2;
%                 macro_dyn.Irf4.pi5 = pi5;
%                 macro_dyn.Irf4.pi10 = pi10;
                macro_dyn.Irf4.c = ct;
                macro_dyn.Irf4.C = Ct;
                end
            end
      end
      
   end
end
