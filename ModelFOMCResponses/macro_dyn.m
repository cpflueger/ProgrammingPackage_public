%% This class contains all parameters and methods related to the macro dynamics of the model
classdef macro_dyn
    properties
      phi,                       % smoothing parameter for consumption 
      g,                         % consumption growth rate
      delta,                     % leverage parameter relating the dividend growth to consumption growth 
      delta_betaport,            % leverage of beta-sorted portfolios
      theta0,                    % persistence of the surplus consumption ratio
      theta1,                    % dependence of the surplus consumption ratio on current output gap
      theta2,                    % dependence of the surplus consumption ratio on lagged output gap
      gamma,                     % parameter controlling utility curvature
      rf,                        % steady state real short-term interest rate at x_t = 0 
      sigma_vec,                 % variances of fundamental shocks, in particular, those are the calibrated
                                 % values for the variances of respectively, u^IS (set to 0), u^PC, u^MP and u^*
                                 % to get empirically relevant values in Table 4 of the paper, for example sigma^PC = 400*sqrt(sigma_vec(2))
      gamma_x,                   % controls the reaction of the interest rate to the output gap
      gamma_pi,                  % governs the response of the interest rate to inflation 
                                 % relative to its target level
      psi,                       % coefficient on (i_t - E(pi_t+1)) in IS equation
      rho_i,                     % governs the influence of past interest rates on current
                                 % interest rates (MP equation)
      rho_pi,                    % (1-rho_pi), rho_pi are, respectively, the coefficient on 
                                 % expected next period output gap and lagged output gap
                                 % in the PC equation
      rhoxm,                     % coefficient on lagged output gap in IS equation
      rhoxp,                     % coefficient on expected next period output gap in IS equation 
      kappa,                     % controls the sensitivity of inflation to the output gap
      P,                         % P, Q are the matrices determining the solution for the dynamics of 
                                 % Y_t, in the form Y_t = P*Y_t-1 + Q*u_t
      Q,                         %
      number_stable,             % number of stable eigenvalues (modulus < 1)
      number_complex,            % number of complex eigenvalues
      eigenvalues,               % eigenvalues used to construct solution                
      eigenvalues_select,        % the three eigenvalues we actually select
      A,                         % Matrix rotating \hat Y_t so that shocks to Z tilde = A* \hat Y_t 
                                 % are independent standard normal and the first element is conditionally
                                 % perfectly correlated with consumption (Appendix A.2.3)
      Ainv,                      % Inverse of A  
      vast,                      % covariance of inflation target shocks with normalized shocks eps
      sigmaperp,                 % standard deviation of orthogonal component of u^* : refer to (61)-(65) in the appendix
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
      solutionNumber,            % Identifies which solution has been selected by ModelPQ82
   end
   
   methods
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% update_params updates the Euler equation parameters that depend on fundamental parameters
       %% and set the variance shocks vector and the    
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function macro_dyn = update_params(macro_dyn)  
           %% Euler equation parameters
           macro_dyn.theta2         = (macro_dyn.phi-1)-macro_dyn.theta1; 
           macro_dyn.rhoxp          = 1/(macro_dyn.phi - macro_dyn.theta1); 
           macro_dyn.rhoxm          = macro_dyn.theta2/(macro_dyn.phi - macro_dyn.theta1); 
           macro_dyn.psi            = 1/(macro_dyn.gamma*(macro_dyn.phi - macro_dyn.theta1));   
           %% Variance shocks vector
           % The volatilities of non-MP shocks (e.g. demand or supply shocks) should not be changed because they may generate wrong results.
           macro_dyn.sigma_vec      = [10^(-9) 10^(-9) macro_dyn.sigma_vec(3) 10^(-9)];  
           assert((macro_dyn.sigma_vec(1) == 10^(-9)) && (macro_dyn.sigma_vec(2) == 10^(-9)) && (macro_dyn.sigma_vec(4) == 10^(-9)),...
                '@update_params: Do not change volatility of non-MP shocks. The asset price solution is not set up to handle this.')
           %% Industry portfolio betas
           macro_dyn.delta_betaport = [1.4813;1.0255;1.0255;0.9132;0.8031;0.7017;0.6802;0.6410;0.5376;0.4796];
       end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% ModelPQ82 Solves for the macro dynamics of the model using the method of generalized eigenvectors and selects an equilibrium 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function macro_dyn = ModelPQ82(macro_dyn, num_set)
            %% Set up matrices F,G,H,M, defining the equation to be solved 
            % 0 = F*E_t(Y_t+1) + G*Y_t + H*Y_t-1 + M*u_t
            F= [macro_dyn.rhoxp, macro_dyn.psi      , 0;
                0              ,(1-macro_dyn.rho_pi) , 0;
                0              , 0                   , 0];
           
            G= [-1                                   , 0                                      ,-macro_dyn.psi;
                macro_dyn.kappa                      , -1                                     ,0             ;
                (1-macro_dyn.rho_i)*macro_dyn.gamma_x, (1-macro_dyn.rho_i)*macro_dyn.gamma_pi ,-1           ];    
            
            H= [macro_dyn.rhoxm, 0               , 0         ;
                0              , macro_dyn.rho_pi, 0               ;
                0              , 0               , macro_dyn.rho_i];
            
            M= [1, 0, 0, 0                ;
                0, 1, 0, -macro_dyn.rho_pi;
                0, 0, 1, -macro_dyn.rho_i];
        
            % G-matrix used to solve for lead-lag coefficients
            G_a=G;
        %% Solve for P and Q
            m=3; % Number of state variables. In our case we have three: 
            % output gap, inflation and nominal interest rate gaps.

            % Selecting a solution is equivalent to picking three out of the 
            % six generalized eigenvalues, hence we impose a series of refinements as in
            % section A.2 of the appendix, recalling that eigenvalues with absolute
            % value <1 give stable solution

            % Construct matrices Xi and Delta as in Uhlig(1999) and compute the
            % generalized eigenvealues and eigenvectors 
            Xi=    [-G_a    , -H       ;
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
            % Count how many stable and complex eigenvalues there are
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
        
            Ps         = zeros(size(select_mat,1),3,3);
            Qs         = zeros(size(select_mat,1),3,4);
            
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

                % Output gap shock loading onto u_t
                QMs(i,:) =[1,0,0]*macro_dyn.Q; 
            end  
        
       % We used the first solution as the best, even though we have 2 or more solutions
       macro_dyn.eigenvalues_select = evec2(select_mat(1,:)); 
       macro_dyn.P          = reshape(Ps(1,:,:),3,3);
       macro_dyn.Q          = reshape(Qs(1,:,:),3,4);
       macro_dyn.QM         = QMs(1,:);
       macro_dyn.solutionNumber = 1;
      end
     
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Compute the state vector Ztilde=A \hat Y
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function macro_dyn = ScaledStateVector(macro_dyn)
         
           % Obtain the variance-covariance matrix of fundamental shocks
           macro_dyn.Sigmau = diag(macro_dyn.sigma_vec);
           
           % QM is loading of consumption shock onto u_t
           macro_dyn.QM     = [1,0,0]*macro_dyn.Q;
                 
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
            % Pick A2 in the null space of A1*Q*Sigmau*Q'
            null1 = null(A1*macro_dyn.Q*macro_dyn.Sigmau*macro_dyn.Q');
            A2 = null1(:,1)';
            % Scale so that A2*u has unit variance
            A2 = A2/sqrt(A2*macro_dyn.Q*macro_dyn.Sigmau*macro_dyn.Q'*A2');
            
            % Pick A3 in null space of A1*Q*Sigmau*Q' and A2*Q*Sigmau*Q'
            null2 = null([A1*macro_dyn.Q*macro_dyn.Sigmau*macro_dyn.Q'; A2*macro_dyn.Q*macro_dyn.Sigmau*macro_dyn.Q']);
            A3=null2';
            A3=A3/sqrt(A3*macro_dyn.Q*macro_dyn.Sigmau*macro_dyn.Q'*A3');
            macro_dyn.A = [A1;A2;A3];
                        
            % Compute A^{-1}
            macro_dyn.Ainv=inv(macro_dyn.A);

            % Compute loading of v_* onto normalized vector of shocks \eps_t. Corresponds
            % to vec_* in the appendix
            macro_dyn.vast=macro_dyn.A*macro_dyn.Q*[0,0,0,1]'*macro_dyn.sigma_vec(4);

            % sigmaperp = standard deviation of orthogonal component of v_* 
            macro_dyn.sigmaperp = sqrt(macro_dyn.sigma_vec(4)-macro_dyn.vast'*macro_dyn.vast);
            
            % Dynamics for scaled state vector \tilde P
            macro_dyn.Ptilde                 = macro_dyn.A*macro_dyn.P/macro_dyn.A;
            
            % Unconditional variance and standard deviation of \tilde Z
            % StdZtilde will be used determine the numerical grid for the value
            % function iteration
            macro_dyn.VarZtilde              = Uncon_var(macro_dyn.Ptilde, diag([1,1,1]));
            macro_dyn.StdZtilde              = sqrt(diag(macro_dyn.VarZtilde));
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% ModelPQ82_dynare Solves for the macro dynamics of the model using Dynare
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function macro_dyn = ModelPQ82_dynare(macro_dyn)
            % Solving the model in dynare
            dynare NK_dynare.mod
            % Loading model results
            load NK_dynare\Output\NK_dynare_results.mat
            
            % Obtaining the stable eigenvalues from the solved model in dynare
            stable_eigenval = find(abs(oo_.dr.eigval)<1);
            
            % Identifying the components of the matrices P and Q from the solved model in dynare
            order_y  = find(oo_.dr.order_var==1);
            order_i  = find(oo_.dr.order_var==2);
            order_pi = find(oo_.dr.order_var==3);

            y_t  = oo_.dr.ghx(order_y,:);
            i_t  = oo_.dr.ghx(order_i,:);
            pi_t = oo_.dr.ghx(order_pi,:);

            e_y_t  = oo_.dr.ghu(order_y,:);
            e_i_t  = oo_.dr.ghu(order_i,:);
            e_pi_t = oo_.dr.ghu(order_pi,:);
           
            % Filling in the macro model results in the macro_dyn class
            macro_dyn.eigenvalues_select    =  oo_.dr.eigval(stable_eigenval); 
            macro_dyn.P                     =  [y_t(2)  y_t(3)  y_t(1);
                                                pi_t(2) pi_t(3) pi_t(1);
                                                i_t(2)  i_t(3)  i_t(1)];
            macro_dyn.Q                     =  [e_y_t(2)  e_y_t(1)  e_y_t(3)  e_y_t(4);
                                                e_pi_t(2) e_pi_t(1) e_pi_t(3) e_pi_t(4);
                                                e_i_t(2)  e_i_t(1)  e_i_t(3)  e_i_t(4)];
            macro_dyn.QM                    =  macro_dyn.Q(1,:);
            
      end

  end

end
