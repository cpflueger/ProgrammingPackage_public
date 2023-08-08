%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% num_set is a class containing all settings for the numerical solution of asset prices
% num_set also calculates the weights used for Gauss-Legendre integration and the grid 
% for the state vector, that we use when solving for asset prices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef num_set
    properties
        Nbonds,         % maximal bond maturity 
        NN,             % maturity of zero coupon bonds/dividend claims 
        Nbins,          % number of bins to make the plot of regression coefficients by bin (10 for deciles etc)
        sfast,          % determines which grid to use 
        N,              % number of gridpoints along each dimension (appendix A.2.4)
        m,              % width of the grid as a multiple of the unconditional standard deviation of Z tilde
        GLpoints,       % number of Gauss Legendre points for integration along dimension 1
        GLpoints2,      % number of Gauss Legendre points for integration along dimension 2
        GLpoints3,      % number of Gauss Legendre points for integration along dimension 3
        GLdomain,       % domain of integration = 8 standard deviations for first step of integration
        GLdomain2,      % domain of integration for dimension 2
        GLdomain3,      % domain of integration for dimension 3
        Nsim,           % Number of simulations (set to >=2)
        T,              % length of simulation
        Tirf,           % length of impulse response
        burn,           % length of burn-in period
        N1,             % number of gridpoints to list in Z
        N2,             % number of gridpoints to list in shat 
        N3,             % number of gridpoints to list in Ztilde 
        sizexm,         % number of gridpoints (stored in xmgrid)for previous period de-trended real rate \hat r_{t-1}
        xmgrid,         % list of gridpoints for previous period de-trended real rate
        splusgrid,      % grid for s
        xGL,            % Gauss-Legendre grid points on dimension 1
        xGL2,           % Gauss-Legendre grid points on dimension 2
        xGL3,           % Gauss-Legendre grid points on dimension 3
        wGL,            % Gauss-Legendre weights for dimension 1
        wGL2,           % Gauss-Legendre weights for dimension 2
        wGL3,           % Gauss-Legendre weights for dimension 3
        prob1,          % scaled weights for Gauss-Legendre on dimension 1
        prob2,          % scaled weights for Gauss-Legendre on dimension 2
        prob3,          % scaled weights for Gauss-Legendre on dimension 3    
        Z,              % list (i.e. matrix in which each row is a point) of (3 dimensional) gridpoints for ztilde
        X,              % list of gridpoints (single coordinate as x is scalar) for x
        Ztilde,         % list of gridpoints in the subspace along dimensions 2 and 3 
                        % i.e. for dimensions 2 and 3 of z tilde 
        shat,           % list of  scalar log deviations of surplus-consumption ratio from steady state
        lambdas,        % values of sensitivity function at each of the points listed in shat
        h               % bond maturity in quarters on RHS of FOMC day yield regressions
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Update some specifications of class num_set
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function num_set = update_num(num_set, macro_dyn)
            % Drop some observations at the beginning of each simulation
            num_set.burn = 100; 
            % Generate the grid for surplus consumption ratio shat
            num_set      = num_set.chooseshatgrid(macro_dyn); 
            % Generates the grid for value function iteration 
            num_set      = num_set.creategrid(macro_dyn);
            % Evaluate lambda at all gridpoints
            num_set                = num_set.senshat(macro_dyn);
            % Generate transition probabilities
            num_set                = num_set.generateprobs;
            % Generate splusgrid 
            num_set                = num_set.generatesplusgrid(macro_dyn);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Fill in some parameters needed to run the num_set class functions 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function num_set = parameters(num_set)
            num_set.N = 2;
            num_set.T = 10000; 
            num_set.Nsim = 1; 
            num_set.sfast = 6; 
            num_set.m = [2 2 2]'; 
            num_set.sizexm = 2; 
            num_set.GLpoints2 = 10; 
            num_set.GLdomain2 = 6; 
            num_set.GLpoints = 40; 
            num_set.GLdomain = 8; 
            num_set.Nbonds = 40;
            num_set.NN = 300; 
            num_set.Nbins = 10; 
            num_set.h = 4; 
        end      
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% creategrid generates the grid for value function iteration 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function num_set = creategrid(num_set, macro_dyn)
            countpoints           = zeros(max(num_set.N));
            countpoints(1,1:num_set.N)    = (1:num_set.N);
            countpoints(2,1:num_set.N)    = (1:num_set.N);
            countpoints(3,1:num_set.N)    = (1:num_set.N);

            % Create grid with dimension length(StdZtilde) X N i.e. 3 X N
            zgrid = -repmat(num_set.m.*macro_dyn.StdZtilde, 1, num_set.N)+ 2*(countpoints-1).*repmat(num_set.m.*macro_dyn.StdZtilde./(num_set.N-1), 1, num_set.N); 

            % Create the matrix Z in which each row is a point in the 3 dimensional space of the
            % Z tilde state variables: essentially listing gridpoints
            [z1, z2, z3] = ndgrid(zgrid(1,:), zgrid(2,:), zgrid(3,:)); 
            z1=reshape(z1,1,num_set.N^3); 
            z2=reshape(z2,1,num_set.N^3);
            z3=reshape(z3,1,num_set.N^3);
            num_set.Z = vertcat(z1,z2,z3)';

            % Create list of the values for the real rate deviations from steady-state used in the grid (recall Ztilde = A*Yhat = A*[x,pihat,ihat]')   
            num_set.X  = ([0,0,1]-[0,1,0]*macro_dyn.P)*macro_dyn.Ainv*num_set.Z';

            % Generate a grid for dimensions 2 and 3 
            % Ztilde lists all possible grid points in the subspace along dimensions 2 and 3
            zgridtilde   = zgrid(2:3,:);
            [z1, z2]     = ndgrid(zgridtilde(1,:), zgridtilde(2,:));
            PN           = num_set.N^2;
            z1           = reshape(z1,1,PN);
            z2           = reshape(z2,1,PN);
            num_set.Ztilde       = vertcat(z1,z2)';
            
            % Create grid for \hat r_{t-1}=r_{t-1}-\bar r
            min_X             = min(num_set.X );
            max_X             = max(num_set.X );
            num_set.xmgrid         = (min_X:(max_X-min_X)/(num_set.sizexm-1):max_X);

            % Generate points and weights for numerical integration
            num_set.GLpoints3              = num_set.GLpoints2;
            num_set.GLdomain3              = num_set.GLdomain2;

            [num_set.xGL, num_set.wGL]     = GaussLegendre(num_set.GLpoints);
            num_set.xGL                    = num_set.xGL*num_set.GLdomain;
            num_set.wGL                    = num_set.wGL*num_set.GLdomain;

            [num_set.xGL2, num_set.wGL2]   = GaussLegendre(num_set.GLpoints2);
            num_set.xGL2                   = num_set.xGL2*num_set.GLdomain2;
            num_set.wGL2                   = num_set.wGL2*num_set.GLdomain2;

            [num_set.xGL3, num_set.wGL3]   = GaussLegendre(num_set.GLpoints3);
            num_set.xGL3                   = num_set.xGL3*num_set.GLdomain3;
            num_set.wGL3                   = num_set.wGL3*num_set.GLdomain3;

            % Number of grid points for Ztilde
            num_set.N1                     = num_set.N^3;
            % Number of grid points for surplus consumption ratio
            num_set.N2                     = max(size(num_set.shat));
            % Number of grid points on subgrid for Ztilde covering only
            % dimensions 2 and 3
            num_set.N3                     = max(size(num_set.Ztilde));
            
            return
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Generate transition probabilities i.e. scaled Gauss-Legendre weights
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function num_set = generateprobs(num_set)
        % Possible transition points are zexp+e_1 xGl
        % This corresponds to Gauss-Legendre integration scaled to give total
        % probability of one
            num_set.prob1   = num_set.wGL .*(1/sqrt(2*pi)).*exp(-0.5*num_set.xGL .^2);
            wGL_            = [num_set.wGL2, num_set.wGL3]; 
            xGL_            = [num_set.xGL2, num_set.xGL3];
            num_set.prob1   = num_set.prob1/sum(num_set.prob1);
            prob_           = wGL_.*(1/sqrt(2*pi)).*exp(-0.5*xGL_.^2);
            sump_           = sum(prob_);
            num_set.prob2   = prob_(:,1)/sump_(1);
            num_set.prob3   = prob_(:,2)/sump_(2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sensitivity function lambda as a function of surplus consumption ratio
        % shat
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function num_set = senshat(num_set,macro_dyn)
            num_set.lambdas    = real((1/macro_dyn.Sbar)*sqrt(1-2*num_set.shat)-1).*(num_set.shat<=macro_dyn.Smax-macro_dyn.sbar);
            return
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate the distribution of surplus consumption ratio shat(t+1) 
        % at current state vector (Z,shat,x^-). This will be used
        % repeatedly to integrate over saht(t+1) conditional on the current
        % state.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function num_set = generatesplusgrid(num_set, macro_dyn)
        
        % We evaluate \tilde u_{1,t} at these 40 points
        % index1: z_t
        % index2: shat_t
        % index3: x_{t-1}
        % index4: e_2 z_{t+1}, e_3 z_{t+1} - we need this to interpolate over dimensions 2 and 3 at time t+1
        % index5: epsilon_{1,t+1}
        
        % Define constants to get how  E_t s_{t+1} loads onto state vector
        phi    =   macro_dyn.phi;
        Ainv   =   macro_dyn.Ainv;
        P      =   macro_dyn.P;
        const2 =  [1,0,0]*( P- phi*eye(3))* Ainv;
        const3 = ([0,0,1]-[0,1,0]* P)* Ainv;

        % Initialize splusgrid 
        InitialValue    = zeros(num_set.N1, num_set.N2, num_set.sizexm, num_set.N3, num_set.GLpoints);
        num_set.splusgrid       = InitialValue;

        % Consumption shock at each of the Gauss-Legendre points for
        % \eps_{1,t+1}
        epsc                    = macro_dyn.sigmac*num_set.xGL;

        % Loop over grid for state vector and compute vector of possible
        % values for s(t+1) for each Gauss-Legendre point of \eps_{1,t+1}
            for i1=1:num_set.N1
                for i4=1:num_set.N3
                    for i5=1:num_set.GLpoints
                        for i2=1:num_set.N2
                            for i3=1:num_set.sizexm                              
                                % Alternative definition, that does not
                                % explicitly rely on x_{t-1}
                                sint  = macro_dyn.theta0*num_set.shat(i2)+((1/macro_dyn.gamma)*const3-const2)*num_set.Z(i1,:)'+num_set.lambdas(i2)*epsc(i5);
                                num_set.splusgrid(i1,i2,i3,i4,i5) =  sint;
                            end
                        end
                    end
                end
            end
        % Ensure that splus and factor2 are real
        num_set.splusgrid   = real(num_set.splusgrid);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate the grid for surplus consumption ratio shat
        % Our baseline is grid 6
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function num_set = chooseshatgrid(num_set, macro_dyn)
            switch(num_set.sfast)
                case {1,2,3,4,5}
                    switch(num_set.sfast)
                        case 1,
                            % Grid 1: Same as original CC grid with 13 points
                            splot  = log(macro_dyn.Smax/13:macro_dyn.Smax/13:macro_dyn.Smax);
                            num_set.shat           = real(splot-macro_dyn.sbar);
                        case 2,
                            % Grid 2: Same as original CC grid with 50 points
                            splot  = log(macro_dyn.Smax/50:macro_dyn.Smax/50:macro_dyn.Smax);
                            num_set.shat           = real(splot-macro_dyn.sbar);
                        case 3,
                            % Grid 2: Same as original CC grid with 100 points
                            splot  = log(macro_dyn.Smax/100:macro_dyn.Smax/100:macro_dyn.Smax);
                            num_set.shat   = real(splot-macro_dyn.sbar);
                        case 4,
                            num_set.shat   = [-296.2995  , -196.2995, -30:10:-2, -1.2995];
                        case 5,
                            num_set.shat   = [-296.2995  , -196.2995, -30:5:-2 , -1.2995];
                    end
                splot= num_set.shat + macro_dyn.sbar;

                case {6,7,8}
                    switch(num_set.sfast)
                        case 6, % This is the case we mostly use
                            % Intermediate size grid of size 50
                            splotu  = log(macro_dyn.Smax/20:macro_dyn.Smax/20:macro_dyn.Smax);
                            % Lower segment as in Wachter
                            splot1  = (-50:(50+splotu)/30:min(splotu)-(50+splotu)/130);
                            splot   = horzcat(splot1, splotu); %concatenate arrays horizontally
                        case 7,
                            splotu  = log(macro_dyn.Smax/20:macro_dyn.Smax/20:macro_dyn.Smax);
                            % Lower segment as in Wachter
                            splot1=(-100:(100+splotu)/100:min(splotu)-(100+splotu)/130);
                            splot  = horzcat(splot1, splotu);
                        case 8,
                            % This is the largest grid, exactly as in Wachter
                            % Upper segment as in Wachter. 
                            splotu  = log(macro_dyn.Smax/101:macro_dyn.Smax/101:macro_dyn.Smax);
                            % Lower segment as in Wachter
                            splot1=(-300:(300+splotu)/900:min(splotu)-(300+splotu)/900);
                            splot  = horzcat(splot1, splotu);
                    end

            % Subtract sbar to obtain grid for surplus consumption ratio relative to steady-state  
            num_set.shat = real(splot-macro_dyn.sbar);
            end
        end
    end
end

