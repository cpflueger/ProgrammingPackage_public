%% This class contains all parameters and methods related to the macro dynamics of the model
classdef macro_dyn
    properties
      phi,                       % smoothing parameter for consumption 
      g,                         % consumption growth rate
      delta,                     % leverage parameter relating the dividend growth to consumption growth 
      theta0,                    % persistence of the surplus consumption ratio
      theta1,                    % dependence of the surplus consumption ratio on current output gap
      theta2,                    % dependence of the surplus consumption ratio on lagged output gap
      gamma,                     % parameter controlling utility curvature
      rf,                        % steady state real short-term interest rate at x_t = 0 
      sigma_vec,                 % variances of orthogonal shocks (u_t)
      sigma_vNew,                % standard deviations of model shocks (v_t)
      corr_vNew,                 % correlations of model shocks (v_t)
      psi,                       % coefficient on (i_t - E(pi_t+1)) in log-linear consumption Euler equation (real rate slope)
      rhoxm,                     % coefficient on lagged output gap in log-linear consumption Euler equation (rho^x in main paper)
      rhoxp,                     % coefficient on expected next period in log-linear consumption Euler equation (f^x in main paper)
      P,                         % P, Q are the matrices determining the solution for the dynamics of 
                                 % Y_t, in the form Y_t = P*Y_t-1 + Q*u_t
      Q,                         %
      number_stable,             % number of stable eigenvalues (modulus < 1)
      number_complex,            % number of complex eigenvalues
      eigenvalues,               % eigenvalues used to construct solution    
      eigenvalues_select,        % the three eigenvalues we actually select
      A,                         % Matrix defining the scaled state vector such that \tilde Z= A* \hat Y_t 
      Ainv,                      % Inverse of A  
      vast,                      % covariance of v_* with normalized shocks eps. This is vec* in the appendix
      sigmaperp,                 % standard deviation of orthogonal component of v_*
      Ptilde,                    % matrix determining dynamics of scaled state vector \tilde Z 
      VarZtilde,                 % Unconditional variance of tilde Z
      StdZtilde,                 % Unconditional standard deviations of \tilde Z (3X1 vector as we have 3 state variables) 
      sigmac,                    % Standard deviation of consumption shocks epsilon_c,t 
      betaq,                     % quarterly pure time discount rate 
      ImpliedParams,             % vector of parameters implied by the solution found: discount rate beta, 
                                 % Euler equation coefficients (rho^x, f^x, \psi), steady-state surplus consumption ratio Sbar, 
                                 % maximum log surplus consumption ratio s^max, exp(s^max), std. of consumption shocks in annualized percent
      Sbar,                      % steady state surplus-consumption ratio  
      sbar,                      % log of Sbar
      smax,                      % value of s (log surplus-consumption ratio) such that for larger s 
                                 % the sensitivity function equals zero (as in Campbell-Cochrane (1999))  
      Smax,                      % exp(smax)
      Sigmau,                    % variance covariance matrix of orthogonal shocks (u_t)
      QM,                        % Consumption loading onto vector of shocks u_t
      J,                         % SMM objective function 
      solutionNumber             % Identifies which solution has been selected by ModelPQ82
   end
   
   methods
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% update_params updates the Euler equation parameters that depend on fundamental parameters
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function macro_dyn = update_params(macro_dyn)
           macro_dyn.rhoxp = 1/(macro_dyn.phi - macro_dyn.theta1);
           macro_dyn.rhoxm = macro_dyn.theta2/(macro_dyn.phi - macro_dyn.theta1);
           macro_dyn.psi    = 1/(macro_dyn.gamma*(macro_dyn.phi - macro_dyn.theta1));
       end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% ModelPQ82 Solves for the macro dynamics of the model using the method of generalized eigenvectors and selects an equilibrium 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function macro_dyn = ModelPQ82(macro_dyn, num_set)
        %define parameters  
        theta1=macro_dyn.theta1;
        theta2=macro_dyn.theta2;
        phi=macro_dyn.phi;
        gamma=macro_dyn.gamma;
        Pinput=macro_dyn.P(2:3,:);

        %Define matrix M that relates model shocks v_t to vector of
        %orthogonal shocks u_t via v_t=M u_t
        m23 = macro_dyn.sigma_vNew(1)*macro_dyn.corr_vNew(3,1)/macro_dyn.sigma_vNew(3);
        m31 = macro_dyn.sigma_vNew(2)*(macro_dyn.corr_vNew(2,1) - macro_dyn.corr_vNew(3,1)*macro_dyn.corr_vNew(3,2))/(macro_dyn.sigma_vNew(1)*(1-macro_dyn.corr_vNew(3,1)^2));
        m33 = macro_dyn.sigma_vNew(2)*macro_dyn.corr_vNew(3,2)/macro_dyn.sigma_vNew(3);
        
        Mtemp = [1,   0, m23;
                 m31, 1, m33; 
                 0,   0, 1];
        
        %variance-covariance matrix of v_t
        Sigma_vNew = corr2cov(macro_dyn.sigma_vNew, macro_dyn.corr_vNew);
        %variance-covariance matrix of u_t
        macro_dyn.sigma_vec = [10^-16,diag(inv(Mtemp)*Sigma_vNew *inv(Mtemp)')'];

       %% Define Blanchard-Kahn problem
       Atilde=[1,1/gamma,0,0;
        0,1,0,0;
        0,0,1,0;
        -(phi-theta1),0,-1/gamma,1];

        Btilde=zeros(4,4);
        Btilde(1,4)=1;
        Btilde(2:3,1:3)=Pinput; 
        Btilde(4,1)=-theta2;

        A=inv(Atilde)*Btilde;

        %eigenvalue decomposition
        [B,J]=eig(A);             
        
        % Sort eigenvalues
        evec=diag(J);
        [~,IX] = sort(abs(evec));  
        evec2  = evec(IX); 
        macro_dyn.eigenvalues = evec2;
        B2     = B(:,IX);  
        
        %Equilibrium selection        
        % detect stable eigenvalues with modulus less than one
        evec_stable = evec2(abs(evec2)<1);
        macro_dyn.number_stable  = size(evec_stable,1);
        % Count how many stable and complex eigenvalues there are
        evec_complex    = evec2(and(real(evec2)~=evec2,abs(evec2)<1)); 
        macro_dyn.number_complex = size(evec_complex,1); 

        % Implement the following with a switch
        % 3 stable evalues > unique solution
        % 4 stable evalues (2 complex) > 2 possible solutions
        % 4 stable evalues (all real) > 4 possible solutions
        
        %select_mat is such that each row corresponds to a selection of
        %three eigenvalues, that lead to a real-valued non-explosive
        %solution

        switch macro_dyn.number_stable
            case 2
                disp('No solutions!')
            case 3
                select_mat = [true,true,true,false];
            case 4
                if macro_dyn.number_complex == 0
                    select_mat = [[true,true,true,false];
                                  [true,true,false,true];
                                  [true,false,true,true];
                                  [false,true,true,true]];
                end
                if macro_dyn.number_complex == 2
                    select_complex = and(real(evec2)~=evec2,abs(evec2)<1);
                    select_complex = select_complex';
                    others = [1,1,1,1] - select_complex;
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

        end
        %compute Ps, Qs, impulse responses and SMM objective function for all
        %real-valued non-explosive equilibria
        Js         = zeros(size(select_mat,1),1);
        Ps         = zeros(size(select_mat,1),3,3);
        Qs         = zeros(size(select_mat,1),3,4);
        QMs        = zeros(size(select_mat,1),4);
        
        %loop over different equilibria
        for i = 1:size(select_mat,1)
            select = select_mat(i,:); 
            
        % Solution for P and Q
        B11=B2(1:3,select);
        J1=diag(evec2(select));
            
            %solution for P
            macro_dyn.P = real(B11*J1*inv(B11)); % This is B^(BK) in the appendix
            % check that B11 is not too close to be not invertible
            if min(abs(eig(B11))) < 1e-6
                warning('Cant construct Q, skipping this combination of eigenvalues')
                Js(i) = 10000;
                continue;
            end    
            %matrix Sigma from BK solution
            Sigma=zeros(3,2);  % This is Sigma^(BK) in the appendix
            Sigma(2,1)=1;
            Sigma(3,2)=1;
           
            if (phi-theta1-macro_dyn.P(1,1)-macro_dyn.P(2,1)/gamma) < 1e-10
                warning('Cant construct Q, skipping this combination of eigenvalues')
                Js(i) = 10000;
                continue;
            end         
                    
            Sigma(1,1)=(macro_dyn.P(1,2)+macro_dyn.P(2,2)/gamma)/(phi-theta1-macro_dyn.P(1,1)-macro_dyn.P(2,1)/gamma);
            Sigma(1,2)=(macro_dyn.P(1,3)+(macro_dyn.P(2,3)-1)/gamma)/(phi-theta1-macro_dyn.P(1,1)-macro_dyn.P(2,1)/gamma);
            %corresponding solution for Q, which designates the loadings
            %onto the orthogonalizes shocks u_t
            macro_dyn.Q      = zeros(3,4);
            
            macro_dyn.Q(:,2:4) = Sigma*Mtemp(1:2,1:3);
            Ps(i,:,:)         = macro_dyn.P;
            Qs(i,:,:)         = macro_dyn.Q;
        end
        
       %We used the first solution as the best, even though we have 2 or more solutions
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
           macro_dyn.QM     = [1,0,0]*macro_dyn.Q-(macro_dyn.phi-macro_dyn.theta1)*[1,0,0,0];
           
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
            % pick A3 in null space of Q
            null2 = null(macro_dyn.Q(:,end-2:end)');
            macro_dyn.A = [A1;A2;null2'];
            
            % Compute A^{-1}
            macro_dyn.Ainv=inv(macro_dyn.A);

            % Compute loading of v_* onto normalized vector of shocks \eps_t. Corresponds
            % to vec_* in the appendix
            macro_dyn.vast=macro_dyn.A*macro_dyn.Q*[0,0,0,1]'*macro_dyn.sigma_vec(4);

            % sigmaperp = standard deviation of orthogonal component of v_* 
            macro_dyn.sigmaperp = sqrt(macro_dyn.sigma_vec(4)-macro_dyn.vast'*macro_dyn.vast);
            
            % dynamics for scaled state vector \tilde P
            macro_dyn.Ptilde                 = macro_dyn.A*macro_dyn.P/macro_dyn.A;
            
            % unconditional variance and standard deviation of \tilde Z
            % StdZtilde will be used determine the numerical grid for the value
            % function iteration
            macro_dyn.VarZtilde              = Uncon_var(macro_dyn.Ptilde, diag([1,1,0]));
            macro_dyn.StdZtilde              = sqrt(diag(macro_dyn.VarZtilde));
      end
  end

end
