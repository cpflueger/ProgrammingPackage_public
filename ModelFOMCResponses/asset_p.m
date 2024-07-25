%% This class contains all parameters and methods related to asset prices and their properties
classdef asset_p
    %% List of properties
    properties
        G,                         % price-consumption ratio of consumption claim, evaluated at each point of the N1 X N2 X N3 grid
        Bn,                        % real bond prices, again evaluated on a grid with same dimensions
                                   % but for 19 different maturities (from 2 up to 20 quarters), hence
                                   % has dimension (Nbonds-1) X N1 X N2 X N3
        Bnom,                      % nominal bond prices, same dimensions as Bn
        G_rn,                      % Those correspond to the three matrices above, computed with risk-neutral pricing kernel
        Bn_rn,                     %
        Bnom_rn,                   %
        risk_neutral_run,          % = 0 or 1, if equal to 1 compute risk neutral asset prices setting gamma = 0
        ImpliedParams,             % vector of parameters implied by the solution found (computed in macro_dyn): discount rate beta,
                                   % Euler equation coefficients (rho^x, f^x, \psi), steady-state surplus consumption ratio Sbar,
                                   % maximum log surplus consumption ratio s^max, exp(s^max), std. of consumption shocks in annualized percent
        rho,                       % constant used in the calculation of real rate news in stock reuturns Campbell-Shiller
        stocks,                    % Structure of moments related to stocks:
                                   %    Equity Premium, Equity volatility,
                                   %    Sharpe ratio, mean price-dividend
                                   %    ratio, volatility price-dividend
                                   %    ratio, AR(1) coefficient
                                   %    price-dividend ratio, coeffcient
                                   %    1-year excess return onto pd,
                                   %    R-squared 1-year excess return onto pd
        nominalBonds,              % Structure of moments related to nominal bonds:
                                   %    Term premium, bond return volatility, Sharpe ratio, mean log yield spread,
                                   %    volatility log yield spread, AR(1) coefficient log yield spread,
                                   %    coefficient 1-year bond excess returns onto lagged log yield spread,
                                   %    R^2 1-year bond excess returns onto lagged log yield spread
                                   %    inflBondVarRatio: ratio of variance of inflation changes divided by 10 year yield bond changes
        realBonds,                 % Structure of moments related to real bonds:
                                   %    Term premium, volatility of real bond returns, Sharpe ratio, mean log yield spread,
                                   %    volatility log yield spread, real bond-stock beta, correlation real bond returns with stock returns
        breakevens,                % Structure of moments related to inflation breakeven:
                                   % Term premium, volatility of breakeven bond returns, Sharpe ratio
        crossAsset,                % Structure of cross-asset moments:
                                   %    correlation nominal bond returns with stock returns
                                   %    Beta of nominal bond return on stock return
                                   %    Correlation quarterly inflation-output gap
                                   %    correlation 5-year average inflation-output gap
                                   %    correlation 5-year average Fed Funds-output gap
                                   %    Cross-asset and macro predictability:
                                   %    Coefficient and R^2 of regression of 1-year excess stock return on output gap,
                                   %    Coefficient and R^2 of regression of 1-year excess bond return on output gap
        
        macroDynamics,             % Structure of moments relating to macro dynamics
                                   %   volatility and AR(1) of nominal short rate changes
                                   %   volatility and AR(1) of inflation changes
                                   %   volatility and AR(1) of log real consumption growth
                                   %   volatility and AR(1) output gap
        
        AP_responses,              %   coefficients of bond and stock returns components (cf, real rate, risk premium) onto monetary policy shocks
        additional_moments,        %   fraction simulation periods with s_t>s^max, std approximation error for 1-period nominal rate
        initialShockVec,           %   this vector is needed for simulating IRFs of macro dynamics and asset prices to various IRFl shocks
        Irf3                       %   impulse responses in response to third structural shock
        Irftemp,                   %   temporary impulse response. This field contains the most recently calculated structural impulse response
        testirf,                   %   equals 1 if running impulse responses in test mode without background noise
        FOMC_run,                  %   equals 1 if simulating FOMC effects
        simulated,                 %   simulated series that can reproduce some important moments 
        simulated_rn,              %   simulated series that can reproduce some important moments when we are in the risk neutral case
        risk_neutral_run_aux,      %   auxiliary dummy variable needed to run simulate_moments correctly
    end
    %% Methods section
    methods
        %% Initialize output asset_p class and save implied parameters to solve for asset prices
        function asset = initial_ap(assetInput,macro_dyn)
            % Initialize output asset_p class
            rng(0)
            asset = asset_p(assetInput.G,assetInput.Bn,assetInput.Bnom,0);
            % Set up for risk premia version (that we use then when evaluate
            % the asset price recursion)
            asset.risk_neutral_run = 0;
            % Save implied parameters
            asset.ImpliedParams    = macro_dyn.ImpliedParams;
            asset.initialShockVec  = assetInput.initialShockVec;
            asset.FOMC_run         = assetInput.FOMC_run;
        end
        %% Solve and simulate risk neutral asset prices
        function asset = risk_neutral_ap(assetInput,macro_dyn,num_set)
            % The conditional is for when you want to get the results for 
            % the neutral risk case
            if assetInput.risk_neutral_run == 1
                asset = assetInput.initial_ap(macro_dyn);
                disp('Computing risk neutral prices')
                assetRN = asset_p(0,0,0,0);
                assetRN.risk_neutral_run = 1;
                tic
                assetRN = assetRN.computeFn21(num_set,macro_dyn);
                toc
                assetRN.G_rn = assetRN.G;
                assetRN.Bn_rn = assetRN.Bn;
                assetRN.Bnom_rn = assetRN.Bnom;
                assetRN.initialShockVec = assetInput.initialShockVec;
                assetRN.FOMC_run = assetInput.FOMC_run;
                % Set a preliminary value for log-linearization constant.
                % it does not enter into any numbers shown in the paper and
                % simply serves to avoid an error message
                assetRN.rho = 0.99;
                % Simulate risk neutral prices and returns
                disp('Simulate risk neutral moments')
                tic
                assetRN = assetRN.SimulateMoments(num_set, macro_dyn);
                toc
                % Risk-neutral price-dividend ratio
                PDRN=assetRN.stocks.meanPDlev;
                % Log-linearization constant to decompose risk-neutral 
                % returns into cash flow news and discount rate news
                asset.rho=1/(1+1/(4*PDRN));
                % Save risk-neutral asset prices in aset
                asset.G_rn = assetRN.G;
                asset.Bn_rn = assetRN.Bn;
                asset.Bnom_rn = assetRN.Bnom;
                asset.simulated_rn = assetRN.simulated;
            else
                asset = assetInput;
                asset.ImpliedParams    = macro_dyn.ImpliedParams;
                asset.initialShockVec  = assetInput.initialShockVec;
                asset.FOMC_run         = assetInput.FOMC_run;
                asset.G       = zeros(num_set.N1,num_set.N2,num_set.sizexm);
                asset.Bn      = zeros(num_set.Nbonds-1,num_set.N1,num_set.N2,num_set.sizexm);
                asset.Bnom    = zeros(num_set.Nbonds-1,num_set.N1,num_set.N2,num_set.sizexm);
                asset.G_rn    = zeros(num_set.N1,num_set.N2,num_set.sizexm);
                asset.Bn_rn   = zeros(num_set.Nbonds-1,num_set.N1,num_set.N2,num_set.sizexm);
                asset.Bnom_rn = zeros(num_set.Nbonds-1,num_set.N1,num_set.N2,num_set.sizexm);
                asset.rho     = 0.9956; 
            end
            asset.risk_neutral_run_aux  = 1;
        end
        %% Constructor method
        function asset = asset_p(G , Bn, Bnom, ImpliedParams)
            % This method initializes the class given some input arguments
            % and fills the remaining ones with default values
            if nargin == 0
                asset.G                  = 0;
                asset.Bn                 = 0;
                asset.Bnom               = 0;
                asset.ImpliedParams      = zeros(1,10);
                asset.additional_moments = zeros(11,1);
                
            else
                asset.G                  = G;
                asset.Bn                 = Bn;
                asset.Bnom               = Bnom;
                asset.ImpliedParams      = ImpliedParams;
                asset.additional_moments = zeros(11,1);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SimulateMoments simulates macroeconomic dynamics and asset prices to
        % obtain moments reported in the paper
        % key OUTPUTS:
        %           asset.macroDynamics = macroeconomic dynamics of model
        %           asset.stocks = moments for stock prices and returns
        %           asset.nominalBonds = moments for nominal bond prices and returns
        %           asset.realBonds = moments for real bond prices and returns
        %           asset.crossAsset = moments for the comovement of bonds and stocks and real and nominal macroeconomic variables
        %           asset.additional moments = fraction of surplus consumption
        %           ratio simulations above smax, std of approximation error for
        %           one-period nominal interest rate
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function asset  = SimulateMoments(asset_pI, num_set, macro_dyn)
            asset_pI.risk_neutral_run   = asset_pI.risk_neutral_run_aux;
            Nsim    = num_set.Nsim;
            T       = num_set.T;
            N       = num_set.N;
            sizexm  = num_set.sizexm;
            shat    = num_set.shat;
            xmgrid  = num_set.xmgrid;
            N2      = num_set.N2;
            burn    = num_set.burn;
            z       = num_set.Z;
            Nbonds  = num_set.Nbonds;
            Nbins   = num_set.Nbins;
            h       = num_set.h;
            gamma     = macro_dyn.gamma;
            delta     = macro_dyn.delta;
            if isempty(macro_dyn.delta_betaport)
                % If the leverage of beta-sorted portfolios is uninitialized,
                % then set it to be the leverage of the market portfolio.
                delta_betaport = macro_dyn.delta;
            else
                delta_betaport = macro_dyn.delta_betaport;
            end
            theta0    = macro_dyn.theta0;
            P         = macro_dyn.P;
            Q         = macro_dyn.Q;
            A         = macro_dyn.A;
            Ainv      = macro_dyn.Ainv;
            Sigmau    = macro_dyn.Sigmau;
            Ptilde    = macro_dyn.Ptilde;
            g         = macro_dyn.g;
            phi       = macro_dyn.phi;
            sigmac    = macro_dyn.sigmac;
            rf        = macro_dyn.rf;
            Sbar      = macro_dyn.Sbar;
            smax      = macro_dyn.smax;
            QM        =[1,0,0]*macro_dyn.Q;
            sigma_vec = macro_dyn.sigma_vec;
            initialShockVec = asset_pI.initialShockVec;
            
            % Pre-allocate constants used for dynamics of shat
            const2 =  [1,0,0]*( P- phi*eye(3))* Ainv;
            const3 = ([0,0,1]-[0,1,0]* P)* Ainv;
            
            % Copy asset prices
            asset = asset_p(asset_pI.G,asset_pI.Bn,asset_pI.Bnom,macro_dyn.ImpliedParams);
            asset.G_rn              = asset_pI.G_rn;
            asset.Bn_rn             = asset_pI.Bn_rn;
            asset.Bnom_rn           = asset_pI.Bnom_rn;
            asset.rho               = asset_pI.rho;
            asset.FOMC_run          = asset_pI.FOMC_run;
            asset.risk_neutral_run  = asset_pI.risk_neutral_run;
            asset.simulated_rn      = asset_pI.simulated_rn;
                       
            % Standard deviation of shock that is released at quarter-end
            FOMC_std = (asset_pI.initialShockVec)/400;
            
            % Points and weights for numerical integration over MP shock
            num_set            = num_set.generateprobs;
            xGL                = num_set.xGL2;
            
            % Define matrices used often
            Z14                = zeros(1,4);
            Z3T                = zeros(3,T);
            ZT1                = zeros(T,1);                     
            
            % Reshape price-consumption ratio
            G5dim              = log(reshape(asset_pI.G, N, N, N,N2,sizexm));
            
            % Reshape nominal bond prices
            % Horizon for RHS yields in real information effect regressions
          
            Pnomplus5dim       = log(reshape(asset_pI.Bnom(end,:,:,:), N ,N ,N, N2, sizexm));
            Pnomminus5dim      = log(reshape(asset_pI.Bnom(end-1,:,:,:), N ,N ,N, N2, sizexm));
            Pnomhq5dim         = log(reshape(asset_pI.Bnom(h-1,:,:,:), N ,N ,N, N2, sizexm));
            Pnom5y5dim         = log(reshape(asset_pI.Bnom(19,:,:,:), N ,N ,N, N2, sizexm));
            
            % Reshape real bond prices
            Pplus5dim          = log(reshape(asset_pI.Bn(end,:,:,:), N, N, N, N2, sizexm));
            Pminus5dim         = log(reshape(asset_pI.Bn(end-1,:,:,:), N, N, N, N2, sizexm));
            P5y5dim            = log(reshape(asset_pI.Bn(19,:,:,:), N ,N ,N, N2, sizexm));
            
           
            % Do the same for risk neutral asset prices. If no rn then just give the same
            
            if asset.risk_neutral_run == 1
                
                G5dim_rn            = log(reshape(asset_pI.G_rn, N, N, N,N2,sizexm));
                
                Pnomplus5dim_rn     = log(reshape(asset_pI.Bnom_rn(end,:,:,:), N ,N ,N, N2, sizexm));
                Pnomminus5dim_rn    = log(reshape(asset_pI.Bnom_rn(end-1,:,:,:), N ,N ,N, N2, sizexm));
                Pnomhq5dim_rn       = log(reshape(asset_pI.Bnom_rn(h-1,:,:,:), N ,N ,N, N2, sizexm));
                Pnom5y5dim_rn       = log(reshape(asset_pI.Bnom_rn(19,:,:,:), N ,N ,N, N2, sizexm));
                
                Pplus5dim_rn        = log(reshape(asset_pI.Bn_rn(end,:,:,:), N, N, N, N2, sizexm));
                Pminus5dim_rn       = log(reshape(asset_pI.Bn_rn(end-1,:,:,:), N, N, N, N2, sizexm));
                P5y5dim_rn          = log(reshape(asset_pI.Bn_rn(19,:,:,:), N ,N ,N, N2, sizexm));
            else
                
                G5dim_rn             = log(reshape(asset_pI.G, N, N, N,N2,sizexm));
                
                Pnomplus5dim_rn      = log(reshape(asset_pI.Bnom(end,:,:,:), N ,N ,N, N2, sizexm));
                Pnomminus5dim_rn     = log(reshape(asset_pI.Bnom(end-1,:,:,:), N ,N ,N, N2, sizexm));
                Pnomhq5dim_rn        = log(reshape(asset_pI.Bnom(h-1,:,:,:), N ,N ,N, N2, sizexm));
                Pnom5y5dim_rn        = log(reshape(asset_pI.Bnom(19,:,:,:), N ,N ,N, N2, sizexm));
                
                Pplus5dim_rn         = log(reshape(asset_pI.Bn(end,:,:,:), N, N, N, N2, sizexm));
                Pminus5dim_rn        = log(reshape(asset_pI.Bn(end-1,:,:,:), N, N, N, N2, sizexm));
                P5y5dim_rn           = log(reshape(asset_pI.Bn(19,:,:,:), N ,N ,N, N2, sizexm));
            end
            
            
            % Reshape grid for n-dimensional interpolation
            z1shape            = reshape(z(:,1), N, N, N);
            z2shape            = reshape(z(:,2), N, N, N);
            z3shape            = reshape(z(:,3), N, N, N);
            z1grid             = reshape(z1shape(:,1,1), N, 1);
            z2grid             = reshape(z2shape(1,:,1), N, 1);
            z3grid             = reshape(z3shape(1,1,:), N, 1);
            
            % Upper and lower bounds of grids
            zgrid_             = [z1grid, z2grid, z3grid];
            zupper             = max(zgrid_)';
            zlower             = min(zgrid_)';
            supper             = max(shat);
            slower             = min(shat);
            xminusupper        = max(xmgrid);
            xminuslower        = min(xmgrid);
            
            % Start Nsim independent simulations
            rng(0) % To avoid simulation noise
            % j counts the independent simulations
            j=1;
            while j<=1    
                %% Draw shocks - we split all shocks into the component before and after the MP date
                % Generate T draws of u_t and compute \epsilon_t from that
                u            = mvnrnd(Z14, diag(sigma_vec), T)';
                % Information released on MP date
                u_FOMC       = mvnrnd(Z14, diag(FOMC_std.^2), T)';
                u_pre        = u-u_FOMC;
                eps_pre      = A*Q*u_pre;
                eps_FOMC     = A*Q*u_FOMC;
                eps          = eps_pre+eps_FOMC;
                uast_pre     = u_pre(4,:);
                uast_FOMC    = u_FOMC(4,:);
                uast         = uast_pre+uast_FOMC;
                
                % Initialize simulated time series, time series with suffix '_pre' denote prices just prior to MP release
                ztildesim      = Z3T;
                shatsim        = ZT1;
                shatsim_pre    = ZT1;
                Yhat           = Z3T;
                Yhat_pre       = Z3T;
                PD             = ZT1;
                PD_pre         = ZT1;
                PD_rn_pre      = ZT1;
                PD_rn          = ZT1;
                Pnomplus       = ZT1;
                Pnomplus_pre   = ZT1;
                Pnomplus_rn_pre= ZT1;
                Pnomplus_rn    = ZT1;
                Pnomminus      = ZT1;
                Pnomminus_rn   = ZT1;
                Pplus          = ZT1;
                Pplus_pre      = ZT1;
                Pplus_rn_pre   = ZT1;
                Pplus_rn       = ZT1;
                Pminus         = ZT1;
                Pminus_rn      = ZT1;
                csim           = ZT1;
                csim_pre       = ZT1;
                piast          = ZT1;
                piast_pre      = ZT1;
                approxErrI     = ZT1;
                Pnomhq         = ZT1;
                Pnomhq_rn      = ZT1;
                Pnomhq_pre     = ZT1;
                Pnomhq_rn_pre  = ZT1;
                P5y            = ZT1;
                P5y_rn         = ZT1;
                P5y_pre        = ZT1;
                P5y_rn_pre     = ZT1;
                Pnom5y         = ZT1;
                Pnom5y_rn      = ZT1;
                Pnom5y_pre     = ZT1;
                Pnom5y_rn_pre  = ZT1;
                
                % Update state vector
                for t=3:T
                    % Dynamics for \tilde Z
                    ztildesim(:,t) = Ptilde*ztildesim(:,t-1)+eps(:,t);
                    
                    % Dynamics for \hat Y
                    Yhat(:,t)      = P*Yhat(:,t-1)+Ainv*eps(:,t);
                    % Dynamics for surplus consumption ratio relative to steady-state
                    shatsim(t)     = ...
                        theta0*shatsim(t-1)+((1/gamma)*const3-const2)*ztildesim(:,t-1)+...
                        senshat(shatsim(t-1), Sbar)*sigmac*eps(1,t);                    
                    % Truncate state variables at upper and lower end of grid, so we can use standard interpolation
                    zinterp        = max(zlower, min(ztildesim(:,t), zupper));
                    sinterp        = max(slower, min(shatsim(t), supper));
                    xminusinterp   = max(xminuslower, min(const3*ztildesim(:,t-1),xminusupper)); 
                    csim(t)        = g+csim(t-1)+(Yhat(1,t)-phi*Yhat(1,t-1));
                    piast(t)       = piast(t-1)+uast(t);
                    
                    % Scaling factors for nominal bonds
                    ePplus    = exp(Nbonds*piast(t));
                    ePminus   = exp((Nbonds-1)*piast(t));
                    eP5y      = exp(20*piast(t));
                    ePhq      = exp(h*piast(t));
                    
                    % Interpolate to obtain price-consumption ratio and real and nominal bond prices with risk premia
                    qNDInterp = eqinterpnND2(z1grid, z2grid, z3grid, shat, xmgrid, ...
                        {G5dim, Pnomplus5dim, Pnomminus5dim, Pplus5dim, Pminus5dim, Pnomhq5dim, P5y5dim, Pnom5y5dim}, ...
                        zinterp(1), zinterp(2), zinterp(3), ...
                        sinterp, xminusinterp, 0);
                    
                    % Scale price-consumption ratio for consumption claim
                    PD(t)          = qNDInterp(1)/4;
                    
                    % n-quarter bond prices
                    Pnomplus(t)    = qNDInterp(2)/ePplus;
                    Pplus(t)       = qNDInterp(4);
                    
                    % n-1-quarter bond prices
                    Pnomminus(t)   = qNDInterp(3)/ePminus;
                    Pminus(t)      = qNDInterp(5);
                    
                    % h quarters bond prices
                    Pnomhq(t) = qNDInterp(6)/ePhq;
                    
                    % 5 year bond prices
                    P5y(t) = qNDInterp(7);
                    Pnom5y(t) = qNDInterp(8)/eP5y;
                    
                    % Interpolate to obtain risk neutral price-consumption ratio and real and nominal bond prices
                    qNDInterp_rn = eqinterpnND2(z1grid, z2grid, z3grid, shat, xmgrid, ...
                        {G5dim_rn, Pnomplus5dim_rn, Pnomminus5dim_rn, Pplus5dim_rn, Pminus5dim_rn, Pnomhq5dim_rn, P5y5dim_rn, Pnom5y5dim_rn}, ...
                        zinterp(1), zinterp(2), zinterp(3), ...
                        sinterp, xminusinterp, 0);
                    
                    %scale risk-neutral price-consumption ratio for consumption claim
                    PD_rn(t)          = qNDInterp_rn(1)/4;
                    
                    %risk-neutral n-quarter bond prices
                    Pnomplus_rn(t)    = qNDInterp_rn(2)/ePplus;
                    Pplus_rn(t)       = qNDInterp_rn(4);
                    
                    %risk-neutral n-1-quarter bond prices
                    Pminus_rn(t)      = qNDInterp_rn(5);
                    Pnomminus_rn(t)   = qNDInterp_rn(3)/ePminus;
                    
                    % h quarters risk neutral bond prices
                    Pnomhq_rn(t) = qNDInterp_rn(6)/ePhq;
                    
                    % 5 year risk neutral bond prices
                    P5y_rn(t) = qNDInterp_rn(7);
                    Pnom5y_rn(t) = qNDInterp_rn(8)/eP5y;
                    
                    % approximation error for 1-period nominal rate
                    approxErrI(t) = .5*([0,1,0]*Q(:,2:4) + [0,0,1])*Sigmau(2:4,2:4)*([0,1,0]*Q(:,2:4) + [0,0,1])'+ ...
                        gamma*(senshat(shatsim(t), Sbar) + 1)*QM(2:4)*Sigmau(2:4,2:4)*([0,1,0]*Q(:,2:4) + [0,0,1])';
                    
                    % point estimate prior to MP shock
                    ztildesim_pre   = Ptilde*ztildesim(:,t-1)+eps_pre(:,t);
                    shatsim_pre(t)     = ...
                        theta0*shatsim(t-1)+((1/gamma)*const3-const2)*ztildesim(:,t-1)+...
                        senshat(shatsim(t-1), Sbar)*sigmac*eps_pre(1,t);
                    Yhat_pre(:,t)   = P*Yhat(:,t-1)+Ainv*eps_pre(:,t);
                    csim_pre(t)     = g+csim(t-1)+(Yhat_pre(1,t)-phi*Yhat(1,t-1));
                    piast_pre(t)    = piast(t-1)+uast_pre(t);
                    % bound state vectors for interpolation
                    shatsim_pre        = max(slower, min(shatsim_pre, supper));
                    ztildesim_pre       = max(zlower, min(ztildesim_pre, zupper));
                    
                    if asset.FOMC_run ==1
                        
                        ePplus_pre  = exp(Nbonds*piast_pre(t));
                        eP5y_pre    = exp(20*piast_pre(t));
                        ePhq_pre    = exp(h*piast_pre(t));
                        
                        % Interpolate to obtain price-consumption ratio and real and nominal bond prices with risk premia
                        qNDInterp = eqinterpnND2(z1grid, z2grid, z3grid, shat, xmgrid, ...
                            {G5dim, Pnomplus5dim, Pnomminus5dim, Pplus5dim, Pminus5dim, Pnomhq5dim, P5y5dim, Pnom5y5dim}, ...
                            ztildesim_pre(1), ztildesim_pre(2), ztildesim_pre(3), ...
                            shatsim_pre(t), xminusinterp, 0);
                        
                        PD_pre(t)           = qNDInterp(1)/4;
                        Pnomplus_pre(t)     = qNDInterp(2)/ePplus_pre;
                        Pplus_pre(t)        = qNDInterp(4);
                        Pnomhq_pre(t)       = qNDInterp(6)/ePhq_pre;
                        P5y_pre(t)          = qNDInterp(7);
                        Pnom5y_pre(t)       = qNDInterp(8)/eP5y_pre;
                                               
                        qNDInterp_rn = eqinterpnND2(z1grid, z2grid, z3grid, shat, xmgrid, ...
                            {G5dim_rn, Pnomplus5dim_rn, Pnomminus5dim_rn, Pplus5dim_rn, Pminus5dim_rn, Pnomhq5dim_rn, P5y5dim_rn, Pnom5y5dim_rn}, ...
                            ztildesim_pre(1), ztildesim_pre(2), ztildesim_pre(3), ...
                            shatsim_pre(t), xminusinterp, 0);
                        
                        PD_rn_pre(t)           = qNDInterp_rn(1)/4;
                        Pnomplus_rn_pre(t)     = qNDInterp_rn(2)/ePplus_pre;
                        Pplus_rn_pre(t)        = qNDInterp_rn(4);
                        Pnomhq_rn_pre(t)       = qNDInterp_rn(6)/ePhq_pre;
                        P5y_rn_pre(t)          = qNDInterp_rn(7);
                        Pnom5y_rn_pre(t)       = qNDInterp_rn(8)/eP5y_pre;
                        
                    end
                end
                % Nbonds/4-year nominal and real bond yields in annualized percent (in logs!)
                y10nom              = -400*log(Pnomplus)/Nbonds;
                y10real             = -400*log(Pplus)/Nbonds;
                y10nom_rn           = -400*log(Pnomplus_rn)/Nbonds;
                y10real_rn          = -400*log(Pplus_rn)/Nbonds;
                % h quarter nominal bonds yields
                yhqnom              = -400*log(Pnomhq)/h;
                yhqnom_rn           = -400*log(Pnomhq_rn)/h;
                % 5 year bond yields
                y5nom               = -400*log(Pnom5y)/20;
                y5nom_rn            = -400*log(Pnom5y_rn)/20;
                y5real              = -400*log(P5y)/20;
                y5real_rn           = -400*log(P5y_rn)/20;
                
                % Real short rate (i minus expected inflation)
                rfr                = 400*([0,0,1]*Yhat-[0,1,0]*P*Yhat+rf);
                
                % Changes around FOMC announcement
                Ret_FOMC               = exp((csim-csim_pre)).*PD./PD_pre;
                reteq_FOMC             = 100*log((1/delta)*Ret_FOMC - ((1-delta)/delta));
                
                % Generate the time series of beta-sorted portfolio returns,
                % matrix dimension T-by-#port.
                delta_temp = repmat(delta_betaport', T, 1);
                Ret_FOMC_temp = repmat(Ret_FOMC, 1, length(delta_betaport));
                ret_betaport_FOMC = 100*log((1./delta_temp).*Ret_FOMC_temp ...
                    - ((1-delta_temp)./delta_temp));
                
                y10nom_pre              = -400*log(Pnomplus_pre)/Nbonds;
                y10nom_FOMC             = y10nom-y10nom_pre;
                
                y10real_pre             = -400*log(Pplus_pre)/Nbonds;
                y10real_FOMC            = y10real-y10real_pre;
                
                y5real_pre              = -400*log(P5y_pre)/20;
                y5real_FOMC             = y5real - y5real_pre;
                
                y5nom_pre               = -400*log(Pnom5y_pre)/20;
                y5nom_FOMC              = y5nom  - y5nom_pre;
                
                yhqnom_pre              = -400*log(Pnomhq_pre)/h;
                yhqnom_FOMC             = yhqnom - yhqnom_pre;
                
                rfr_nom                 = 400*([0,0,1]*Yhat+piast'+rf);
                rfr_nom_pre             = 400*([0,0,1]*Yhat_pre+piast_pre'+rf);
                rfr_nom_FOMC            = rfr_nom-rfr_nom_pre;
                
                rfr_Pnom                 = 400*([0,0,1]*P*Yhat+piast'+rf);
                rfr_Pnom_pre             = 400*([0,0,1]*P*Yhat_pre+piast_pre'+rf);
                rfr_Pnom_FOMC            = rfr_Pnom-rfr_Pnom_pre;
                                              
                % Risk-neutral changes around FOMC announcement
                y10nom_rn_pre           = -400*log(Pnomplus_rn_pre)/Nbonds;
                y10nom_rn_FOMC          = y10nom_rn - y10nom_rn_pre;
                
                y10real_rn_pre          = -400*log(Pplus_rn_pre)/Nbonds;
                y10real_rn_FOMC         = y10real_rn - y10real_rn_pre;
                
                y5real_rn_pre           = -400*log(P5y_rn_pre)/20;
                y5real_rn_FOMC          = y5real_rn - y5real_rn_pre;
                
                y5nom_rn_pre           = -400*log(Pnom5y_rn_pre)/20;
                y5nom_rn_FOMC          = y5nom_rn - y5nom_rn_pre;
                
                yhqnom_rn_pre           = -400*log(Pnomhq_rn_pre)/h;
                yhqnom_rn_FOMC          = yhqnom_rn-yhqnom_rn_pre;
                
                Ret_rn_FOMC            = exp((csim-csim_pre)).*PD_rn./PD_rn_pre;
                reteq_rn_FOMC          = 100*log((1/delta)*Ret_rn_FOMC- ((1-delta)/delta));
                
                % Risk-neutral beta-sorted portfolio returns,
                % matrix dimension T-by-#port.
                delta_temp = repmat(delta_betaport', T, 1);
                Ret_rn_FOMC_temp = repmat(Ret_rn_FOMC, 1, length(delta_betaport));
                ret_rn_betaport_FOMC = 100*log((1./delta_temp).*Ret_rn_FOMC_temp ...
                    - ((1-delta_temp)./delta_temp));
                                              
                % Levered dividends: Compute D^{delta}_{t+1}
                dDelta = PD(2:end).*exp(csim(2:end)) + exp(csim(2:end)) - (1-delta).*PD(1:end-1).*exp(csim(1:end-1)).*exp(rfr(1:end-1)./400)' ...
                    - delta.*PD(2:end).*exp(csim(2:end));
                
                % Take the 64-quarter moving average of D^{delta}_{t+1}
                dDelta_bar = conv(dDelta,ones(1,64),'valid')/64;
                
                if burn>200
                    dDelta_bar = 0.5*(dDelta_bar(1:end-64)+dDelta_bar(65:end));
                end

                % Drop burn period observations and other observations to ensure that returns and prices have
                % the same length
                PD                 = PD(burn-1+2:end);
                PD_rn              = PD_rn(burn-1+2:end);
                Yhat               = Yhat(:,burn-1+2:end);
                shatsim            = shatsim(burn-1+2:end);
                Pnomplus           = Pnomplus(burn-1+2:end);
                Pnomminus          = Pnomminus(burn-1+2:end);
                Pplus              = Pplus(burn-1+2:end);
                Pminus             = Pminus(burn-1+2:end);
                Pnomplus_rn        = Pnomplus_rn(burn-1+2:end);
                Pnomminus_rn       = Pnomminus_rn(burn-1+2:end);
                Pplus_rn           = Pplus_rn(burn-1+2:end);
                Pminus_rn          = Pminus_rn(burn-1+2:end);
                csim               = csim(burn-1+2:end);
                piast              = piast(burn-1+2:end)';
                y10nom             = y10nom(burn-1+2:end);
                y10real            = y10real(burn-1+2:end);
                rfr                = rfr(burn-1+2:end);
                eps                = eps(:,burn-1+3:end);
                approxErrI         = approxErrI(burn-1+2:end); 
                burn1              = burn+2; 

                if burn>200
                    dDelta_bar         = dDelta_bar((burn-129)+2:end); 
                else
                    dDelta_bar         = dDelta_bar((burn-65)+2:end); 
                end
                
                % Price-dividend ratio of levered equity at time t divides by smoothed dividends
                PDlev              = delta.*PD.*exp(csim)./dDelta_bar;
                
                % 1-quarter nominal yield
                rfr_nom            = 400*([0,0,1]*Yhat+piast+rf);
                
                % Level return on consumption claim
                Ret                = exp((csim(2:end)-csim(1:end-1))).*(1+4*PD(2:end))./(4*PD(1:end-1));
                Ret_rn             = exp((csim(2:end)-csim(1:end-1))).*(1+4*PD_rn(2:end))./(4*PD_rn(1:end-1));
                
                % Log return on levered equity
                reteq              = log((1/delta)*Ret - ((1-delta)/delta)*exp(rfr(1:end-1)'/400));
                reteq_rn           = log((1/delta)*Ret_rn - ((1-delta)/delta)*exp(rfr(1:end-1)'/400));
                delta_temp = repmat(delta_betaport', length(Ret), 1);
                Ret_temp = repmat(Ret, 1, length(delta_betaport));
                rfr_temp = repmat(rfr(1:end-1)'/400, 1, length(delta_betaport));
                ret_betaport = log((1./delta_temp).*Ret_temp - ((1-delta_temp)./delta_temp).*exp(rfr_temp)); 
                
                % In percentage units in excess of the riskfree rate
                reteq              = 100*reteq - rfr(1:end-1)'/4;
                reteq_rn           = 100*reteq_rn - rfr(1:end-1)'/4;
                ret_betaport = 100*ret_betaport - rfr(1:end-1)'/4; %XXXFL
                
                % Log excess nominal bond returns in percent but not annualized
                retnom             = 100*log(Pnomminus(2:end))-100*log(Pnomplus(1:end-1))-rfr_nom(1:end-1)'/4;
                retnom_rn          = 100*log(Pnomminus_rn(2:end))-100*log(Pnomplus_rn(1:end-1))-rfr_nom(1:end-1)'/4;
                
                % Log excess real bond returns in percent but not annualized
                retreal            = 100*log(Pminus(2:end))-100*log(Pplus(1:end-1))-rfr(1:end-1)'/4;
                retreal_rn         = 100*log(Pminus_rn(2:end))-100*log(Pplus_rn(1:end-1))-rfr(1:end-1)'/4;
                
                % Breakeven: log returns on nominal in excess of log returns on real bonds
                breakeven          = retnom - retreal;

                % Nominal and real log yield spreads
                spreadNom          	= y10nom'-rfr_nom;
                spreadReal          = y10real'-rfr;

                % 1-year log equity excess returns in natural units
                ret1yr             = conv(reteq,ones(1,4),'valid')/100;
                
                % Levered price-dividend ratio
                pdlev = log(PDlev);
                
                % Compute cash-flow news and real rate news vectors to compute equity real rate news analytically according to Campbell and Ammer
                rho                = asset_pI.rho;
                rhoIPinv           = inv(eye(3)-rho*P);
                Gammaeq_rr         = -rho*([0,0,1]-[0,1,0]*P)*rhoIPinv*Ainv;
                
                % Cash-flow news of stock and bond returns
                reteq_cf           = reteq_rn-(100*Gammaeq_rr*eps)';
                retnom_cf          = retnom_rn-retreal_rn;
                
                % Risk premium excess returns of stock and bond returns
                reteq_rp           = reteq-reteq_rn;
                retnom_rp          = retnom-retnom_rn;
                                            
%% Regress FOMC asset returns onto nominal rate changes on FOMC dates

                reteq_FOMC     = reteq_FOMC(burn+2:end);
                ret_betaport_FOMC = ret_betaport_FOMC(burn+2:end,:); 
                ret_rn_betaport_FOMC = ret_rn_betaport_FOMC(burn+2:end,:);
                y10nom_FOMC    = y10nom_FOMC(burn+2:end);
                y10real_FOMC   = y10real_FOMC(burn+2:end);
                y5real_FOMC    = y5real_FOMC(burn+2:end);
                y5nom_FOMC     = y5nom_FOMC(burn+2:end);
                yhqnom_FOMC    = yhqnom_FOMC(burn+2:end);
                rfr_nom_FOMC   = rfr_nom_FOMC(burn+2:end);
                rfr_Pnom_FOMC  = rfr_Pnom_FOMC(burn+2:end);
                
                reteq_rn_FOMC   = reteq_rn_FOMC(burn+2:end);
                y10nom_rn_FOMC  = y10nom_rn_FOMC(burn+2:end);
                y10real_rn_FOMC = y10real_rn_FOMC(burn+2:end);
                y5real_rn_FOMC  = y5real_rn_FOMC(burn+2:end);
                y5nom_rn_FOMC   = y5nom_rn_FOMC(burn+2:end);
                yhqnom_rn_FOMC  = yhqnom_rn_FOMC(burn+2:end);
                                
                reteq_rp_FOMC = reteq_FOMC - reteq_rn_FOMC;
                
                bernankeKuttnerTemp1 = regress(reteq_FOMC,[ones(size(rfr_nom_FOMC)); rfr_nom_FOMC]');
                bernankeKuttnerTemp2 = regress(reteq_rn_FOMC,[ones(size(rfr_nom_FOMC)); rfr_nom_FOMC]');
                bernankeKuttner = [bernankeKuttnerTemp1(2),...
                    bernankeKuttnerTemp2(2)];
                
                bernankeKuttnerTemp1 = regress(reteq,[ones(size(rfr_nom_FOMC)); rfr_nom_FOMC]');
                bernankeKuttnerTemp2 = regress(reteq_rn,[ones(size(rfr_nom_FOMC)); rfr_nom_FOMC]');
                bernankeKuttner_month = [bernankeKuttnerTemp1(2),...
                    bernankeKuttnerTemp2(2)];
                

                bernankeKuttnerTemp1 = regress(reteq_FOMC,[rfr_nom_FOMC; rfr_nom_FOMC.*(rfr_nom_FOMC>0); (rfr_nom_FOMC>0); ones(size(rfr_nom_FOMC))]');
                bernankeKuttnerTemp2 = regress(reteq_rn_FOMC,[ rfr_nom_FOMC; rfr_nom_FOMC.*(rfr_nom_FOMC>0); (rfr_nom_FOMC>0); ones(size(rfr_nom_FOMC))]');
                bernankeKuttner_dummy = [bernankeKuttnerTemp1,bernankeKuttnerTemp2];
                               
                % Regressions with portfolio returns on LHS.
                bernankeKuttnerTemp = cell2mat(cellfun(@(y) regress(y, ...
                    [ones(size(y10nom_FOMC)), rfr_nom_FOMC']), ...
                    num2cell(ret_betaport_FOMC, 1), ...
                    'UniformOutput', false));
                bernankeKuttner_betaport = bernankeKuttnerTemp(2,:);
                               
                bernankeKuttnerTemp = cell2mat(cellfun(@(y) regress(y, ...
                    [ones(size(y10nom_FOMC)), rfr_nom_FOMC']), ...
                    num2cell(ret_rn_betaport_FOMC, 1), ...
                    'UniformOutput', false));
                bernankeKuttner_rn_betaport = bernankeKuttnerTemp(2,:);
                  
                % Run regression by deciles
                regCoeffsFF = zeros(Nbins,3);
                
                decilePoints = quantile(shatsim, linspace(0, 1, Nbins + 1));
                for l = 2:(Nbins + 1)
                    relevantPoints = (shatsim < decilePoints(l)) & (shatsim >= decilePoints(l-1));
                    relevantPoints = relevantPoints(2:(T-burn));
                    regTemp = regress(reteq_FOMC(relevantPoints),[ones(sum(relevantPoints),1), rfr_nom_FOMC(relevantPoints)']);
                    regCoeffsFF(l-1,1) = regTemp(2);
                    regTemp = regress(reteq_rn_FOMC(relevantPoints),[ones(sum(relevantPoints),1), rfr_nom_FOMC(relevantPoints)']);
                    regCoeffsFF(l-1,2) = regTemp(2);
                    regTemp = regress(reteq_rp_FOMC(relevantPoints),[ones(sum(relevantPoints),1), rfr_nom_FOMC(relevantPoints)']);
                    regCoeffsFF(l-1,3) = regTemp(2);
                end
                
                regCoeffsQuantilesFF = regCoeffsFF;
                
                % Run regression above and below median
                regCoeffsMedian = zeros(2,3);
                
                    relevantPoints = (shatsim <= decilePoints(6));
                    relevantPoints = relevantPoints(2:(T-burn));
                    regTemp = regress(reteq_FOMC(relevantPoints),[ones(sum(relevantPoints),1), rfr_nom_FOMC(relevantPoints)']);
                    regCoeffsMedian(1,1) = regTemp(2);
                    regTemp = regress(reteq_rn_FOMC(relevantPoints),[ones(sum(relevantPoints),1), rfr_nom_FOMC(relevantPoints)']);
                    regCoeffsMedian(1,2) = regTemp(2);
                    regCoeffsMedian(1,3) = regCoeffsMedian(1,1)-regCoeffsMedian(1,2);

                    relevantPoints = (shatsim > decilePoints(6));
                    relevantPoints = relevantPoints(2:(T-burn));
                    regTemp = regress(reteq_FOMC(relevantPoints),[ones(sum(relevantPoints),1), rfr_nom_FOMC(relevantPoints)']);
                    regCoeffsMedian(2,1) = regTemp(2);
                    regTemp = regress(reteq_rn_FOMC(relevantPoints),[ones(sum(relevantPoints),1), rfr_nom_FOMC(relevantPoints)']);
                    regCoeffsMedian(2,2) = regTemp(2);
                    regCoeffsMedian(2,3) = regCoeffsMedian(2,1)-regCoeffsMedian(2,2);

                % Run regression with median interaction
                    lowS = (shatsim <= decilePoints(6));
                    lowS=  lowS(2:(T-burn));
                    lowS=  lowS';
                    highS=(lowS==0);

                bernankeKuttnerTemp1 = regress(reteq_FOMC,[rfr_nom_FOMC; rfr_nom_FOMC.*lowS; lowS; ones(size(rfr_nom_FOMC))]');
                bernankeKuttnerTemp2 = regress(reteq_rn_FOMC,[rfr_nom_FOMC; rfr_nom_FOMC.*lowS; lowS; ones(size(rfr_nom_FOMC))]');
                bernankeKuttner_interact = [bernankeKuttnerTemp1,...
                    bernankeKuttnerTemp2];
                
                bernankeKuttnerTemp1 = regress(reteq_FOMC(lowS),[rfr_nom_FOMC(lowS); rfr_nom_FOMC(lowS)]');
                bernankeKuttnerTemp2 = regress(reteq_rn_FOMC(lowS),[rfr_nom_FOMC(lowS); rfr_nom_FOMC(lowS)]');
                bernankeKuttner_lowS = [bernankeKuttnerTemp1,...
                bernankeKuttnerTemp2];
             
                bernankeKuttnerTemp1 = regress(reteq_FOMC(highS),[rfr_nom_FOMC(highS); rfr_nom_FOMC(highS)]');
                bernankeKuttnerTemp2 = regress(reteq_rn_FOMC(highS),[rfr_nom_FOMC(highS); rfr_nom_FOMC(highS)]');
                bernankeKuttner_highS = [bernankeKuttnerTemp1,...
                bernankeKuttnerTemp2];

               % Portfolios and median interaction
               bernankeKuttnerTemp = cell2mat(cellfun(@(y) regress(y, ...
                    [rfr_nom_FOMC', rfr_nom_FOMC'.*lowS', lowS',ones(size(y10nom_FOMC))]), ...
                    num2cell(ret_betaport_FOMC, 1), ...
                    'UniformOutput', false));
                bernankeKuttner_betaport_interact = bernankeKuttnerTemp;
                               
                bernankeKuttnerTemp = cell2mat(cellfun(@(y) regress(y, ...
                    [ones(size(y10nom_FOMC)), rfr_nom_FOMC']), ...
                    num2cell(ret_rn_betaport_FOMC, 1), ...
                    'UniformOutput', false));
                bernankeKuttner_rn_betaport = bernankeKuttnerTemp(2,:);
                              
                realInfoEffectTemp1 = regress(y5real_FOMC,[ones(size(yhqnom_FOMC)), rfr_nom_FOMC']);
                realInfoEffectTemp2 = regress(y5real_rn_FOMC,[ones(size(yhqnom_FOMC)), rfr_nom_FOMC']);
                realInfoEffectTemp3 = regress(y5real_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                realInfoEffectTemp4 = regress(y5real_rn_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                realInfoEffectTemp5 = regress(y10real_FOMC,[ones(size(yhqnom_FOMC)), rfr_nom_FOMC']);
                realInfoEffectTemp6 = regress(y10real_rn_FOMC,[ones(size(yhqnom_FOMC)), rfr_nom_FOMC']);
                realInfoEffectTemp7 = regress(y10real_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                realInfoEffectTemp8 = regress(y10real_rn_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                
                realyieldsMP = [realInfoEffectTemp1(2),...
                    realInfoEffectTemp5(2),...
                    realInfoEffectTemp3(2),...
                    realInfoEffectTemp7(2);...
                    realInfoEffectTemp2(2),...
                    realInfoEffectTemp6(2),...
                    realInfoEffectTemp4(2),...
                    realInfoEffectTemp8(2)];
                
                nomInfoEffectTemp1 = regress(y5nom_FOMC,[ones(size(yhqnom_FOMC)), rfr_nom_FOMC']);
                nomInfoEffectTemp2 = regress(y5nom_rn_FOMC,[ones(size(yhqnom_FOMC)), rfr_nom_FOMC']);
                nomInfoEffectTemp3 = regress(y5nom_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                nomInfoEffectTemp4 = regress(y5nom_rn_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                nomInfoEffectTemp5 = regress(y10nom_FOMC,[ones(size(yhqnom_FOMC)), rfr_nom_FOMC']);
                nomInfoEffectTemp6 = regress(y10nom_rn_FOMC,[ones(size(yhqnom_FOMC)), rfr_nom_FOMC']);
                nomInfoEffectTemp7 = regress(y10nom_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                nomInfoEffectTemp8 = regress(y10nom_rn_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                
               nomyieldsMP = [nomInfoEffectTemp1(2),...
                    nomInfoEffectTemp5(2),...
                    nomInfoEffectTemp3(2),...
                    nomInfoEffectTemp7(2);...
                    nomInfoEffectTemp2(2),...
                    nomInfoEffectTemp6(2),...
                    nomInfoEffectTemp4(2),...
                    nomInfoEffectTemp8(2)];
                
                realInfoEffectTemp1 = regress(y5real_FOMC,[ones(size(yhqnom_FOMC)), rfr_Pnom_FOMC']);
                realInfoEffectTemp2 = regress(y5real_rn_FOMC,[ones(size(yhqnom_FOMC)), rfr_Pnom_FOMC']);
                realInfoEffectTemp3 = regress(y5real_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                realInfoEffectTemp4 = regress(y5real_rn_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                realInfoEffectTemp5 = regress(y10real_FOMC,[ones(size(yhqnom_FOMC)), rfr_Pnom_FOMC']);
                realInfoEffectTemp6 = regress(y10real_rn_FOMC,[ones(size(yhqnom_FOMC)), rfr_Pnom_FOMC']);
                realInfoEffectTemp7 = regress(y10real_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                realInfoEffectTemp8 = regress(y10real_rn_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                
                realyieldsMP_forward = [realInfoEffectTemp1(2),...
                    realInfoEffectTemp5(2),...
                    realInfoEffectTemp3(2),...
                    realInfoEffectTemp7(2);...
                    realInfoEffectTemp2(2),...
                    realInfoEffectTemp6(2),...
                    realInfoEffectTemp4(2),...
                    realInfoEffectTemp8(2)];
                
                nomInfoEffectTemp1 = regress(y5nom_FOMC,[ones(size(yhqnom_FOMC)), rfr_Pnom_FOMC']);
                nomInfoEffectTemp2 = regress(y5nom_rn_FOMC,[ones(size(yhqnom_FOMC)), rfr_Pnom_FOMC']);
                nomInfoEffectTemp3 = regress(y5nom_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                nomInfoEffectTemp4 = regress(y5nom_rn_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                nomInfoEffectTemp5 = regress(y10nom_FOMC,[ones(size(yhqnom_FOMC)), rfr_Pnom_FOMC']);
                nomInfoEffectTemp6 = regress(y10nom_rn_FOMC,[ones(size(yhqnom_FOMC)), rfr_Pnom_FOMC']);
                nomInfoEffectTemp7 = regress(y10nom_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                nomInfoEffectTemp8 = regress(y10nom_rn_FOMC,[ones(size(yhqnom_rn_FOMC)), yhqnom_FOMC]);
                
               nomyieldsMP_forward = [nomInfoEffectTemp1(2),...
                    nomInfoEffectTemp5(2),...
                    nomInfoEffectTemp3(2),...
                    nomInfoEffectTemp7(2);...
                    nomInfoEffectTemp2(2),...
                    nomInfoEffectTemp6(2),...
                    nomInfoEffectTemp4(2),...
                    nomInfoEffectTemp8(2)];
                                                
                %% Stock moments
                
                % Equity risk premium in annualized units
                EqPremium       = 4*(mean(reteq) + .5*std(reteq)^2/100);
                
                % Std. log equity excess returns
                Stdeq           = std(reteq)*2;
                
                % Exp(mean(log pd))
                mean_pdlev      = exp(mean(pdlev));
                
                % Standard deviation of log dp
                std_dp          = std(pdlev);
                
                % Autocorrelation of dp
                dp_corr            = corrcoef(pdlev(2:end), pdlev(1:end-1));
                rho_dp          = dp_corr(1,2);
                
                % Predictability with pd 1 quarter, 1 year and 5 year regressions of returrn on price-dividend ratios
                [re1_coef,~,~,~,R2_re1] = regress(ret1yr, [ones(size(pdlev(1:end-4))), pdlev(1:end-4)]);
                re1                  = re1_coef(2);
                re1_r2               = R2_re1(1);
                
                %% Nominal bond moments
                % Nominal term premium
                BondPremium     = 4*(mean(retnom) + .5*std(retnom)^2/100);
                % Std of bond returns
                Stdnom          = std(retnom)*2;
                
                % Mean and std of log yield spread
                TermSlope       = mean(spreadNom);
                TermSlopeStd    = std(spreadNom);
                
                % Autocorrelation of log yield spread
                AR_slope5temp      = corrcoef(spreadNom(2:end), spreadNom(1:end-1));
                AR_slope5       = AR_slope5temp(1,2);
                
                % Regress returns onto lagged log yield spread
                ret1yrNom          = retnom(4:end)+retnom(3:end-1)+retnom(2:end-2)+retnom(1:end-3);
                ret1yrNom          = ret1yrNom/100;
                % Multiply returns by 100 to match units in empirical exercise
                [ys1_coef,~,~,~,R2_ys1] = regress(100*ret1yrNom, [ones(size(spreadNom(1:end-4)')), spreadNom(1:end-4)']);
                ys1             = ys1_coef(2);
                ys1_r2          = R2_ys1(1);
                
                %bond returns are predicted by PD ratio
                [ypd1_coef,~,~,~,R2_ypd1] = regress(100*ret1yrNom, [ones(size(spreadNom(1:end-4)')), pdlev(1:end-4)]);
                ypd1             = ypd1_coef(2);
                ypd1_r2          = R2_ypd1(1);
                
                [ypd2_coef,~,~,~,R2_ypd2] = regress(100*(ret1yrNom(1:end-4)+ret1yrNom(5:end)), [ones(size(spreadNom(1:end-8)')), pdlev(1:end-8)]);
                ypd2             = ypd2_coef(2);
                ypd2_r2          = R2_ypd2(1);
                
                %% Cross-asset
                % Bond-stock return correlations
                bondstock_corr_temp     = corrcoef(retnom, reteq);
                tipsstock_corr_temp     = corrcoef(retreal, reteq);
                correlations       = [bondstock_corr_temp(1,2), tipsstock_corr_temp(1,2)];
                
                % Nominal bond beta
                beta_temp           = regress(retnom, [ones(T-burn1+1,1), reteq]);
                beta_nom         = beta_temp(2);
                
                % 1-year excess stock return on output gap
                [coeffStockGap_temp,~,~,~,R2StockGap_temp]         = regress(ret1yr, [ones(size(Yhat(1,1:end-4)')), Yhat(1,1:end-4)']);
                coeffStockGap   = coeffStockGap_temp(2);
                R2StockGap      = R2StockGap_temp(1);
                
                % 1-year excess bond return on output gap
                [coeffBondGap_temp,~,~,~,R2BondGap_temp]         = regress(ret1yrNom, [ones(size(Yhat(1,1:end-4)')), Yhat(1,1:end-4)']);
                coeffBondGap    = coeffBondGap_temp(2);
                R2BondGap       = R2BondGap_temp(1);
                
                % Betas of beta-sorted portfolios.
                beta_temp = cell2mat(cellfun(@(y) regress(y, [ones(length(reteq),1), reteq]), ...
                    num2cell(ret_betaport, 1), ...
                    'UniformOutput', false));
                beta_port = beta_temp(2, :)';
                
                %% Real bonds moments
                % Term premium
                RealBondPremium = 4*(mean(retreal) + .5*std(retreal)^2/100);
                % Std returns
                Stdreal         = std(retreal)*2;
                
                % Mean and std. log yield spread
                TermSlopeReal   = mean(spreadReal);
                TermSlopeRealStd= std(spreadReal);
                % Real bond beta
                beta_temp          = regress(retreal, [ones(T-burn1+1,1), reteq]);
                beta_real       = beta_temp(2);
                
                %% Breakeven moments
                BreakevenTermPremium      =  4*(mean(breakeven) + .5*std(breakeven)^2/100);
                StdBreakeven              =  std(breakeven)*2;
                SharpeRatioBreakeven      =  BreakevenTermPremium(j)/StdBreakeven(j);
                
                beta_temp                    = regress(breakeven, [ones(T-burn1+1,1), reteq]);
                beta_breakeven            = beta_temp(2);
                breakevenSim(1:(T-burn1+2))  = y10nom - y10real;
                %% Macro dynamics
                % Std and AR(1) of changes in nominal 1-quarter yield
                rfrNomStd       = std(rfr_nom(5:end)-rfr_nom(1:end-4));
                
                % Std inflation changes
                piChanges            = 4*(100*Yhat(2,1:end-4)+100*piast(1:end-4) - (100*Yhat(2,5:end)+100*piast(5:end)));
                piChangeVol       = std(piChanges);
                         
                % Std Log consumption growth
                consGrowthVol     = 100*std(csim(5:end)-csim(1:end-4));
                
                % Std Output gap
                xVol              = 100*std(Yhat(1,5:end)-Yhat(1,1:end-4));
                
                % Slope coefficient consumption growth output gap growth
                slopetemp            = regress((csim(5:end)-csim(1:end-4)), [ones(size(Yhat(1,5:end))); Yhat(1,5:end)-Yhat(1,1:end-4)]');
                slope_cx             = slopetemp(2,1);

                % Correlation(Delta rfr, Delta consumption growth)
                corrtemp                 = corrcoef(csim(5:end)'-csim(1:end-4)', rfr_nom(5:end)-rfr_nom(1:end-4));
                corr_c_rfr            = corrtemp(1,2);
                
                % Fraction s_t>s^max
                fracmax         = sum(shatsim+log(Sbar)>smax)/T;
                % Std of approximation error for 1-period nominal rate in annualized percent
                stdApproxErrI     = std(approxErrI)*400;
                
                % Skip simulation run if PDlev turns negative
                if min(PDlev)<0
                    j=j+1;
                    continue
                end
                % Increase loop counter
                j=j+1;
            end
                      
            %% Save output
            asset.additional_moments     = [mean(fracmax), stdApproxErrI];
            
            %% Equities
            stocks.equityPremium     = EqPremium;
            stocks.vol               = mean(Stdeq);
            stocks.sharpeRatio       = EqPremium/mean(Stdeq);
            stocks.meanPDlev         = mean_pdlev;
            stocks.stdDP             = std_dp;
            stocks.rhoDP             = rho_dp;
            stocks.coeffRegRetOnPD1y = re1;
            stocks.R2RegRetOnPD1y    = re1_r2;
            
            asset.stocks             = stocks;
            
            %% Nominal Bonds
            nominalBonds.termPremium               = BondPremium;
            nominalBonds.vol                       = mean(Stdnom);
            nominalBonds.sharpeRatio               = BondPremium/mean(Stdnom);
            nominalBonds.meanLogYieldSpread        = TermSlope;
            nominalBonds.volLogYieldSpread         = TermSlopeStd;
            nominalBonds.persistenceLogYieldSpread = AR_slope5;
            nominalBonds.betaNom                   = beta_nom;
            nominalBonds.coeffRegRetOnYS1y         = ys1;
            nominalBonds.R2RegRetOnYS1y            = ys1_r2;
            nominalBonds.coeffRegRetOnPD1y         = ypd1;
            nominalBonds.R2RegRetOnPD1y            = ypd1_r2;
            nominalBonds.coeffRegRetOnPD2y         = ypd2;
            nominalBonds.R2RegRetOnPD2y            = ypd2_r2;
            
            asset.nominalBonds                = nominalBonds;
            
            %% Real Bonds
            realBonds.termPremium          = RealBondPremium;
            realBonds.vol                  = mean(Stdreal);
            realBonds.sharpeRatio          = RealBondPremium/mean(Stdreal);
            realBonds.meanLogYieldSpread   = TermSlopeReal;
            realBonds.volLogYieldSpread    = TermSlopeRealStd;
            realBonds.betaRealStock        = beta_real;
            realBonds.corrRealStock        = correlations(2);
            
            asset.realBonds = realBonds;
            
            %% Brekaven moments
            breakevens.termPremium          = BreakevenTermPremium;
            breakevens.vol                  = StdBreakeven;
            breakevens.sharpeRatio          = SharpeRatioBreakeven;
            breakevens.stockBeta            = beta_breakeven;
            breakevens.simulation           = breakevenSim;
            asset.breakevens = breakevens;
            
            %% Cross-Asset moments
            crossAsset.corrNomStock     = correlations(1);
            crossAsset.betaNom          = beta_nom;
            crossAsset.coeffStockGap    = coeffStockGap;
            crossAsset.R2StockGap       = R2StockGap;
            crossAsset.coeffBondGap     = coeffBondGap;
            crossAsset.R2BondGap        = R2BondGap;
            crossAsset.beta_port        = beta_port; %XXXFL
            
            asset.crossAsset = crossAsset;
            
            %% Macro Dynamics moments
            macroDynamics.iChangeVol     = rfrNomStd;
            macroDynamics.corr_c_rfr     = corr_c_rfr;    
            macroDynamics.piChangeVol    = mean(piChangeVol);
            macroDynamics.consGrowthVol  = mean(consGrowthVol);
            macroDynamics.xVol           = mean(xVol);
            macroDynamics.slope_cx       = mean(slope_cx);
            
            asset.macroDynamics = macroDynamics;
            
            %% Asset returns onto monetary policy shocks
            AP_responses.std_rfr                    = std(rfr_nom_FOMC);    
            %AP_responses.bernankeKuttner            = bernankeKuttner;
            AP_responses.bernankeKuttner            = [bernankeKuttner bernankeKuttner(1)-bernankeKuttner(2)];
            AP_responses.bernankeKuttner_month      = bernankeKuttner_month;
            %AP_responses.bernankeKuttner_dummy      = bernankeKuttner_dummy;
            AP_responses.bernankeKuttner_dummy      = [bernankeKuttner_dummy bernankeKuttner_dummy(:,1)-bernankeKuttner_dummy(:,2)];
            AP_responses.bernankeKuttner_interact   = bernankeKuttner_interact;
            AP_responses.bernankeKuttner_lowS       = bernankeKuttner_lowS;
            AP_responses.bernankeKuttner_highS      = bernankeKuttner_highS;

            AP_responses.bernankeKuttner_betaport               = bernankeKuttner_betaport;
            AP_responses.bernankeKuttner_betaport_interact      = bernankeKuttner_betaport_interact;
            AP_responses.bernankeKuttner_rn_betaport            = bernankeKuttner_rn_betaport;
            AP_responses.betaport_scaled                        = bernankeKuttner_betaport/bernankeKuttner(1);
            AP_responses.betaport_rn_scaled                     = bernankeKuttner_rn_betaport/bernankeKuttner(2);
            
            AP_responses.regCoeffsQuantilesFF        = regCoeffsQuantilesFF;
            AP_responses.regCoeffsMedian             = regCoeffsMedian;
            
            AP_responses.realyieldsMP   = realyieldsMP;
            AP_responses.nomyieldsMP    = nomyieldsMP;
            
            AP_responses.realyieldsMP_forward   = realyieldsMP_forward;
            AP_responses.nomyieldsMP_forward    = nomyieldsMP_forward;
                       
            asset.AP_responses = AP_responses;  
            %% Simulated series

            simulated.pdlev         = pdlev;
            simulated.Yhat          = Yhat;
            simulated.retnom        = retnom;
            simulated.reteq         = reteq;
            simulated.y10nom        = y10nom;
            simulated.rfr_nom       = rfr_nom;  
            simulated.retreal       = retreal;
            simulated.piast         = piast;
            simulated.csim          = csim;
            simulated.reteq_FOMC    = reteq_FOMC;
            simulated.rfr_nom_FOMC  = rfr_nom_FOMC;
            simulated.reteq_rn_FOMC = reteq_rn_FOMC;

            asset.simulated     = simulated;
            return
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % computeFn21 implements the value function iteration
        % key OUTPUTS:
        %   asset.G = price-consumption ratio for claim to all future
        %   consumption on grid points [N1^3 x N2 x sizexm]
        %   asset.Bn = real bond prices for maturities 2 to Nbonds on grid
        %   points
        %   [Nbonds-1 x N1^3 N2 x sizexm]
        %   asset.Bnom = nominal bond prices scaled by exp(-n \pi_{*t}) for
        %   maturities 2 to Nbonds on grid points [Nbonds-1 x N1^3 N2 x
        %   sizexm]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function asset = computeFn21(asset_pI, num_set, macro_dyn)
            %% Initialize output asset_p class and extract needed elements from input classes
            asset = asset_pI;
            risk_neutral_run = asset.risk_neutral_run;
            Nbonds = num_set.Nbonds;
            NN = num_set.NN;
            N = num_set.N;
            sizexm = num_set.sizexm;
            Z = num_set.Z;
            Ztilde = num_set.Ztilde;
            X  = num_set.X ;
            shat = num_set.shat;
            xmgrid = num_set.xmgrid;
            splusgrid = num_set.splusgrid;
            lambdas = num_set.lambdas;
            prob1 = num_set.prob1;
            prob2 = num_set.prob2;
            prob3 = num_set.prob3;
            N1 = num_set.N1;
            N2 = num_set.N2;
            N3 = num_set.N3;
            xGL = num_set.xGL;
            xGL2 = num_set.xGL2;
            xGL3 = num_set.xGL3;
            GLpoints = num_set.GLpoints;
            GLpoints2 = num_set.GLpoints2;
            GLpoints3 = num_set.GLpoints3;
            gamma = macro_dyn.gamma;
            theta0 = macro_dyn.theta0;
            Q = macro_dyn.Q;
            Ainv = macro_dyn.Ainv;
            Sigmau = macro_dyn.Sigmau;
            Ptilde = macro_dyn.Ptilde;
            g = macro_dyn.g;
            phi = macro_dyn.phi;
            sigmac = macro_dyn.sigmac;
            vast = macro_dyn.vast;
            rf = macro_dyn.rf;
            QM = macro_dyn.QM;
            P = macro_dyn.P;
            sigmaperp = macro_dyn.sigmaperp;
            delta = macro_dyn.delta;
            
            % Initialize matrices for log asset prices bnew=log(Bn),
            % bnnew=log(Bnom), fnew=log(Fn)=log price-consumption ratio of zero-coupon
            % consumption claim
            bnew  = zeros(N1,N2, sizexm);
            bnnew = zeros( N1, N2, sizexm);
            fnew  = zeros( N1, N2, sizexm);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Evaluate two-period bond prices and one-period
            % consumption-claim using analytic expressions
            
            % Implement analytic expressions for one-period zero-coupon
            % consumption claim
            
            % Asset prices with risk premia
            if risk_neutral_run == 0
                % Define temporary constants outside the loop to calculate them only once
                const1 =  g;
                const2 =  [1,0,0]*( P- phi*eye(3))* Ainv;
                const3 = ([0,0,1]-[0,1,0]* P)* Ainv;
                % i1 loops over grid points for \tilde Z
                for i1=1: N1
                    % i2 loops over grid points for \hat s_t
                    for i2=1: N2
                        % i3 loops over grid points for x_{t-1}
                        for i3=1: sizexm
                            fnew(i1,i2,i3)= const1 +const2*Z(i1,:)' - rf - const3*Z(i1,:)' - 0.5*gamma*(1- theta0)*(1-2* shat(i2)) + 0.5*(gamma* lambdas(i2)+ (gamma- 1))^2*sigmac^2;
                        end
                    end
                end
                
                % Risk-neutral asset prices
            else
                % Define temporary constants for risk-neutral asset prices outside the loop to calculate them only once
                const1 =  g;
                const2 =  [1,0,0]*( P- phi*eye(3))* Ainv;
                const3 = ([0,0,1]-[0,1,0]* P)* Ainv;
                % i1 loops over grid points for \tilde Z
                for i1=1: N1
                    % i2 loops over grid points for \hat s_t
                    for i2=1: N2
                        % i3 loops over grid points for x_{t-1}
                        for i3=1: sizexm
                            fnew(i1,i2,i3)= const1 + const2*Z(i1,:)' - rf - const3*Z(i1,:)' + 0.5*sigmac^2;
                        end
                    end
                end
            end
            
            % Make sure values are real
            fnew    = real(fnew);
                       
            % Implement analytic expressions for two-period real and
            % nominal bonds
            
            % Define temporary constants outside the loop to calculate them only once
            const1 = -([0,0,1]-[0,1,0]* P)*(eye(3)+ P)* Ainv;
            const2 = -[0,0,1]*(eye(3)+ P)* Ainv;
            vr2     = ([0,0,1]-[0,1,0]* P)* Q;
            vbn2    = [0,1,1]*Q+2*[0,0,0,1];
            
            % Asset prices with risk premia
            if risk_neutral_run == 0
                % i1 loops over grid points for \tilde Z
                for i1 = 1:N1
                    % i2 loops over grid points for \hat s_t
                    for i2 = 1:N2
                        % i3 loops over grid points for x_{t-1}
                        for i3=1: sizexm
                            bnew(i1,i2,i3)  = const1* Z(i1,:)'+0.5*vr2 * Sigmau * vr2' + gamma * (1+lambdas(i2)) * QM * Sigmau * vr2' - 2 * rf;
                            bnnew(i1,i2,i3) = const2* Z(i1,:)'+0.5*vbn2* Sigmau*vbn2'+gamma*(1+lambdas(i2))*QM*Sigmau*vbn2'-2* rf;
                        end
                    end
                end
                
                % Risk-neutral asset prices
            else
                % i1 loops over grid points for \tilde Z
                for i1 = 1:N1
                    % i2 loops over grid points for \hat s_t
                    for i2 = 1:N2
                        % i3 loops over grid points for x_{t-1}
                        for i3 = 1:sizexm
                            bnew(i1,i2,i3)  = const1* Z(i1,:)' + 0.5 * vr2 * Sigmau * vr2' - 2 * rf;
                            bnnew(i1,i2,i3) = const2* Z(i1,:)' + 0.5 * vbn2* Sigmau * vbn2' - 2 * rf;
                        end
                    end
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Prepare inputs to speed up the evaluation of the
            % numerical expectation outside the loop that evaluates the
            % recursion
            
            % Distribution of \tilde Z_{t+1} conditional on \tilde Z_t
            % initialize distributions for three dimensions of \tilde
            % Z_{t+1}
            zp               = zeros(N1, GLpoints);
            zp2              = zeros(N1, GLpoints2);
            zp3              = zeros(N1, GLpoints3);
            % i1 loops over grid for \tilde Z_t
            for i1=1: N1
                % E \tilde Z_{t+1} conditional on Z_t
                ztmp         =  Ptilde* Z(i1,:)';
                % Distribution of first dimension of \tilde Z_{t+1} conditional on \tilde Z_t
                zp(i1,:)     = [1,0,0]*ztmp+ xGL;
                % Distribution of second dimension of \tilde Z_{t+1} conditional on \tilde Z_t
                zp2(i1,:)    = [0,1,0]*ztmp+ xGL2;
                % Distribution of first dimension of \tilde Z_{t+1} conditional on \tilde Z_t.
                % Notice that zp3 is constant across its third dimension, because the shock 
                % \epsilon_{3,t+1} has zero variance by construction.
                zp3(i1,:)    = [0,0,1]*ztmp+ xGL3;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Create map between Z and Ztilde.
            
            % Ai(i1,i4)=1 if and only if the last two dimensions of Z(i1,:)
            % agree with Ztilde(i4,:)
            Ai                = zeros( N1,  N3);
            % i4 loops over Ztilde
            for i4=1: N3
                ztildetemp    =  Ztilde(i4,:);
                % Pick out elements of  Z that correspond to ztildetemp
                Ai(:,i4)      = (...
                    ( Z(:,2)==ztildetemp(:,1)).*...
                    ( Z(:,3)==ztildetemp(:,2))>0);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Prep grid and probabilities so we can take expectations over
            % dimensions 2 and 3 of \epsilon_{t+1} in one step by taking the
            % product of probabilities times realized values along grid for
            % \tilde Z_{2,t+1} and \tilde Z_{3,t+1}
            
            % z1grid lists grid points for \tilde Z_{1,t+1} [N x 1]
            X1           = reshape( Z(:,1), N,  N,  N);
            z1grid       = reshape(X1(:,1,1), N,1);
            
            % X2 lists all gridpoints for \tilde Z_{2,t+1} in [N x N]
            % format
            % X3 lists all grid points for \tilde Z_{3,t+1} in [N x N]
            % format
            % X2 varies along rows and X3 varies along columns
            X2           = reshape( Ztilde(:,1), N,  N);
            X3           = reshape( Ztilde(:,2), N,  N);
            
            % Multiply p2 and p3 to obtain probability weights for the
            % two-dimensional vector (\epsilon_{2,t+1}, \epsilon_{3,t+1}). \tilde
            % Z_{2,t+1} varies along rows of probtilde and \tilde Z_{3,t+1}
            % varies along columns of probtilde
            [p2, p3]     = meshgrid( prob2,  prob3);
            probtilde    = p2.*p3;
            
            % Split up vec*\epsilon_{t+1}=vec*(1)\epsilon_{1,t+1}
            % +vec*(2)\epsilon_{2,t+1}+vec*(3)\epsilon_{3,t+1}
            % uast1, uast2, uast3 save the distribution of the three
            % components. This preallocation again serves to improve speed
            % within the recursion
            uast1        =  vast(1)* xGL';
            uast2        =  vast(2)* xGL2';
            uast3        =  vast(3)* xGL3';
            
            % ua_mesh=vec*(2)\epsilon_{2,t+1}+vec*(3)\epsilon_{3,t+1}
            % \epsilon_{2,t+1} varies along rows and \epsilon_{3,t+1}
            % varies along columns to match the probabilities probtilde
            [ua2, ua3]   = meshgrid(uast2, uast3);
            ua_mesh      = ua2 + ua3;
            
            % Pre-allocate consumption shock distribution for speed
            % consumption shock driven by \epsilon_{1,t+1}
            ve           = [1,0,0]* Ainv;
            ve1          = ve(1)* xGL';
            
            % Pre-allocate inflation shock distribution for speed
            vpi          = [0,1,0]* Ainv;
            vpi1         = vpi(1)* xGL';
            % Inflation shocks driven by \epsilon_{2,t+1} and
            % \epsilon_{3,t+1}. vpi_mesh varies with \epsilon_{2,t+1} along 
            % rows and with \epsilon_{3,t+1} along columns
            vpi2         = vpi(2)* xGL2';
            vpi3         = vpi(3)* xGL3';
            [vpia2, vpia3] = meshgrid(vpi2, vpi3);
            vpi_mesh       = vpia2+vpia3;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Prep interpolation along grid for r^- to improve speed
            % within the main recursion
            
            % For each value in X(i1) weightlow(i1) is such that
            % X(i1)=weightlow(i1)*xmgrid(1)+(1-weightlow(i1)*xmgrid(2)
            mapxlow=zeros( N1,1);
            mapxhigh=zeros( N,1);
            weightlow=zeros( N1,1);
            for i1=1: N1
                findbound       = find(X (i1)>= xmgrid);
                mapxlow(i1)     = max(findbound);
                mapxhigh(i1)    = min(mapxlow(i1)+1, sizexm);
                if mapxlow(i1)~=mapxhigh(i1)
                    weightlow(i1)=1-( X (i1)- xmgrid(mapxlow(i1)))/( xmgrid(mapxhigh(i1))- xmgrid(mapxlow(i1)));
                else
                    weightlow(i1)=1;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Pre-allocate constants outside loop
            const1 =  [1,0,0]*( P- phi*eye(3))* Ainv;
            const2 =  g;
            const3 = ([0,0,1]-[0,1,0]* P)* Ainv;
            const4 = -0.5* gamma*(1- theta0);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Initialize matrices used to extract asset prices along
            % sub-dimensions of the grid within the loop to help with speed and
            % memory
            InitialSheet     = zeros( N2, sizexm);
            InitCondit       = zeros( N3,1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Initialize Price-dividend ratio
            asset_pI.G            = exp(fnew);
            
            % Initialize bond prices
            asset_pI.Bn           = ones(Nbonds-1, N1, N2, sizexm);
            asset_pI.Bnom         = ones(Nbonds-1, N1, N2, sizexm);
            
            % Save 2-period bond prices
            asset_pI.Bn(1,:,:,:)      = exp(bnew);
            asset_pI.Bnom(1,:,:,:)    = exp(bnnew);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Evaluate asset pricing recursion
            % n loops over the maturity of zero coupon bonds/consumption claims
            % The maturity of the zero-coupon consumption claim equals n and the maturity of the
            % zero coupon bonds equals n+1
            
            % A note on notation within the recursion: Variables ending on _old
            % refer to (n-1)-period zero coupon consumption claims and n-period zero
            % coupon bonds. Variables ending on _plus refer to distributions
            % at time t+1.
            for n=2: NN
                % Log of n-1 period price-dividend ratio, real bond price and nominal bond price
                f_old      = fnew;
                bn_old      = bnew;
                bnom_old     = bnnew;
                
                % i1 loops over the grid of current states \tilde Z_t
                for i1=1: N1
                    % Pre-allocate price-dividend ratios conditional on current state \tilde Z_t
                    sheet    = InitialSheet;
                    sheetB   = InitialSheet;
                    sheetBn  = InitialSheet;
                    
                    % i2 loops over grid of current surplus consumption ratio gap \hat s_t
                    for i2=1: N2
                        
                        % i3 loops over lagged detrended real rate hat r_{t-1}
                        for i3=1: sizexm
                            % Distribution of time t+1 log surplus
                            % consumption ratio
                            s_plus           = reshape( splusgrid(i1,i2,i3,1,:)  ,  1,  GLpoints);
                            
                            % Pre-allocate conditional expectations
                            % Fconditional equals the argument in the
                            % recursive expression for F_n,t after taking
                            % expectations over the shock \epsilon_{1,t+1}.
                            % Bconditional and Bnconditional are similarly
                            % the arguments for real and nominal bond
                            % recursions after taking expectations over \epsilon_{1,t+1}.
                            Fconditional    = InitCondit;
                            Bconditional    = InitCondit;
                            Bnconditional   = InitCondit;
                            
                            % i4 loops over points in Ztilde, i.e. the grid
                            % for dimensions 2 and 3
                            for i4=1: N3
                                % Find points in Z corresponding to
                                % Ztilde(i4,:)
                                a     = Ai(:,i4);
                                
                                % Interpolate log-linearly to evaluate f_{n-1,t+1}(\tilde Z_{t+1},s_{t+1},x_t)
                                % Interpolate f_{n-1,t+1}(\tilde Z_{t+1},s_{t+1},x_t) over x_t
                                fmesh_old = reshape(f_old(a>0,:,mapxlow(i1)),  N,  N2)'*weightlow(i1)+reshape(f_old(a>0,:,mapxhigh(i1)),  N,  N2)'*(1-weightlow(i1));
                                % Distribution for \tilde Z_{1,t+1}
                                z1pi  = real(zp(i1,:));
                                
                                % Prepare interpolation/extrapolation of f_{n-1,t+1}(\tilde
                                %Z_{1,t+1},s_{t+1},x_t) over \tilde
                                %Z_{1,t+1} and s_{t+1}
                                M           = size(z1grid,1);
                                [~,XIPOS1]  = histc(min(z1pi,max(X)),z1grid);
                                XIPOS1      = max(XIPOS1,1);
                                XIPOS1      = min(XIPOS1,M-1);
                                T1          =(z1pi'-z1grid(XIPOS1))./(z1grid(XIPOS1+1)-z1grid(XIPOS1));
                                
                                M           = size(shat',1);
                                [~,XIPOS2]  = histc(min(s_plus,max(shat)),shat');
                                XIPOS2      = max(XIPOS2,1);
                                XIPOS2      = min(XIPOS2,M-1);
                                T2          =(s_plus-shat(XIPOS2))./(shat(XIPOS2+1)-shat(XIPOS2));
                                
                                % Interpolate f_{n-1,t+1}(\tilde
                                % Z_{1,t+1},s_{t+1},r_t) over \tilde
                                % Z_{1,t+1}, i.e. first dimension of
                                % \tilde Z_{t+1}
                                FS1         = (1-T1).*diag(fmesh_old(XIPOS2, XIPOS1))+T1.*diag(fmesh_old(XIPOS2, XIPOS1+1));
                                FS2         = (1-T1).*diag(fmesh_old(XIPOS2+1, XIPOS1))+T1.*diag(fmesh_old(XIPOS2+1, XIPOS1+1));
                                
                                % Interpolate f_{n-1,t+1}(\tilde
                                % Z_{1,t+1},s_{t+1},r_t) over s_{t+1}
                                fs_plus          = (1-T2').*FS1+T2'.*FS2;
                                

                                    % This is the recursion equation for
                                    % consumption claims as detailed in
                                    % Appendix section D.1.3
                                    if risk_neutral_run == 0
                                        fs_plus = fs_plus' + const2  + const1 * Z(i1,:)' - rf - const3 * Z(i1,:)' + const4 * (1 - 2* shat(i2)) - (gamma * (1+ lambdas(i2))- 1) * ve1;
                                    else
                                        fs_plus = fs_plus' + const2 +  const1 * Z(i1,:)' - rf - const3 * Z(i1,:)' + ve1;
                                    end
                                    Fconditional(i4)      = real(exp(fs_plus)* prob1);
                                
                                % For n<= maximum bond maturity do the same
                                % calucation to obtain Bconditional and
                                % Bnconditional
                                if n < Nbonds
                                    % Interpolate log-linearly to evaluate
                                    % b_{n-1,t+1}(\tilde
                                    % Z_{t+1},s_{t+1},x_t) and b^\$_{n-1,t+1}(\tilde
                                    % Z_{t+1},s_{t+1},x_t) over x_t
                                    bmesh_old             = reshape(bn_old(a>0,:,mapxlow(i1)),  N,  N2)'*weightlow(i1)+reshape(bn_old(a>0,:,mapxhigh(i1)),  N,  N2)'*(1-weightlow(i1));
                                    Bnmesh_old            = reshape(bnom_old(a>0,:,mapxlow(i1)), N,  N2)'*weightlow(i1)+reshape(bnom_old(a>0,:,mapxhigh(i1)),  N,  N2)'*(1-weightlow(i1));
                                    
                                    % Interpolate b_{n-1,t+1}(\tilde
                                    % Z_{1,t+1},s_{t+1},x_t) over \tilde
                                    % Z_{1,t+1} and s_{t+1}
                                    FS1         = (1-T1).*diag(bmesh_old(XIPOS2, XIPOS1))+T1.*diag(bmesh_old(XIPOS2, XIPOS1+1));
                                    FS2         = (1-T1).*diag(bmesh_old(XIPOS2+1, XIPOS1))+T1.*diag(bmesh_old(XIPOS2+1, XIPOS1+1));
                                    bs_plus          = (1-T2').*FS1+T2'.*FS2;
                                    
                                    % Interpolate b^\$_{n-1,t+1}(\tilde
                                    % Z_{1,t+1},s_{t+1},x_t) over \tilde
                                    % Z_{1,t+1} and s_{t+1}
                                    FS1         = (1-T1).*diag(Bnmesh_old(XIPOS2, XIPOS1))+T1.*diag(Bnmesh_old(XIPOS2, XIPOS1+1));
                                    FS2         = (1-T1).*diag(Bnmesh_old(XIPOS2+1, XIPOS1))+T1.*diag(Bnmesh_old(XIPOS2+1, XIPOS1+1));
                                    bns_plus         = (1-T2').*FS1+T2'.*FS2;
                                    
                                    
                                    % These are the recursion equations for
                                    % real and nominal bonds as detailed in
                                    % Appendix section D.1.3
                                    if risk_neutral_run == 0
                                        bs_plus = bs_plus' - rf - const3 * Z(i1,:)' -0.5 * gamma * (1- theta0) * (1 - 2 * shat(i2)) - gamma * (1 + lambdas(i2)) * ve1;
                                        bns_plus = bns_plus' - rf - [0,0,1] * Ainv * Z(i1,:)' - 0.5 * gamma * (1- theta0) * (1 - 2*shat(i2)) - ...
                                            gamma * (1+ lambdas(i2)) * ve1 - vpi1 - (n+1) * uast1 + 0.5 * (n+1)^2 * sigmaperp^2;
                                    else
                                        bs_plus = bs_plus' - rf - const3 * Z(i1,:)';
                                        bns_plus = bns_plus' - rf -[0,0,1] * Ainv * Z(i1,:)' - vpi1 - (n+1) * uast1 + 0.5 * (n+1)^2 * sigmaperp^2;
                                    end
                                    Bconditional(i4)  = real(exp(bs_plus)* prob1);
                                    Bnconditional(i4) = real(exp(bns_plus)* prob1);
                                end
                            end
                            
                            % Now take expectations of Fconditional,
                            % Bconditional, Bnconditional over the last two
                            % dimensions of \tilde Z_{t+1}
                            
                            % Distributions of \tilde Z_{2,t+1} and \tilde
                            % Z_{3,t+1} conditional on \tilde Z_t
                            zp2i                      = real(zp2(i1,:)');
                            zp3i                      = real(zp3(i1,:)');
                            
                            % Take expectation of Fconditional over \tilde
                                sheet(i2,i3) = expectedinterpolated(Fconditional, zp2i, zp3i, probtilde,  N, X2, X3);
                            
                            if n< Nbonds
                                % Take expectation over Bconditional over \tilde
                                % Z_{2,t+1} and \tilde Z_{3,t+1}
                                sheetB(i2,i3) = expectedinterpolated(Bconditional, zp2i, zp3i, probtilde,  N, X2, X3);
                                
                                % Before taking expectations over Bnconditional
                                % multiply probabilities by the shock to bond
                                % cash flows due to inflation that is perfectly
                                % correlated with \epsilon_{2,t+1} and
                                % \epsilon_{3,t+1}
                                probnom               = probtilde.*exp(-vpi_mesh-(n+1)*ua_mesh);
                                sheetBn(i2,i3)          = expectedinterpolated(Bnconditional, zp2i, zp3i, probnom,  N, X2, X3);
                            end
                        end
                    end
                    fnew(i1,:,:)   = min(max(log(sheet),-300),300);
                    bnew(i1,:,:)   = min(max(log(sheetB),-300),300);
                    bnnew(i1,:,:)  = min(max(log(sheetBn),-300),300);
                end
                
                % Add up zero-coupon consumption claims
                asset_pI.G                  = asset_pI.G+exp(fnew);
                if n< Nbonds
                    % Update real and nominal bond prices
                    asset_pI.Bn(n,:,:,:)    = exp(bnew);
                    asset_pI.Bnom(n,:,:,:)  = exp(bnnew);
                end
            end
            %% Update asset_p output class
                asset.G    = asset_pI.G;
                asset.Bn   = asset_pI.Bn;
                asset.Bnom = asset_pI.Bnom;
            return
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % expectedinterpolated calculates expectations of a function Fconditional(\tilde Z_{2,t+1}, \tilde Z_{3,t+1})
            % It is an intermediate function used in computeFn21
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % INPUTS: Fconditional = known function values at grid points
            %         zp2i = distribution for \tilde Z_{2,t+1}
            %         zp3i = distribution for \tilde Z_{3,t+1}
            %         probtilde = probability that \tilde Z_{2,t+1} = zp2i and
            %         \tilde Z_{3,t+1} = zp3i
            %         N = number of grid points along each dimension \tilde
            %         Z_{2,t+1} and \tilde Z_{3,t+1}
            %         X2 = grid for \tilde Z_{2,t+1} at which function values
            %         are known
            %         X3 = grid for \tilde Z_{3,t+1} at which function values
            %         are known
            
            function Expected  = expectedinterpolated(Fconditional, zp2i, zp3i, probtilde, N, X2, X3)
                
                % Bound Fconditional away from zero to avoid error messages
                Fconditional                 = max(Fconditional, 10^(-320));
                % Take log - we interpolate log-linearly
                fcond                        = log(reshape(Fconditional, [N,N]));
                z2grid                       = X2(:,1,1);
                z3grid                       = X3(1,:,1)';
                
                % Bilinear interpolation over zp3i holding zp2i constant
                % This interpolation is inlined to increase performance
                VORGSIZE = size(fcond);
                Ds       = VORGSIZE(2:end);
                PDs      = prod(Ds);
                APDs     = (1:PDs)';
                fcond        = reshape(fcond,[VORGSIZE(1), PDs]);
                XExt     = {cast(z2grid(:),'double'),APDs};
                F        = griddedInterpolant(XExt,fcond,'linear');
                finterp2       = cast(reshape(F({cast(zp2i(:),class(XExt{1})),APDs}),[length(zp2i), Ds]),superiorfloat(z2grid,fcond,zp2i));
                tmp                          = finterp2';
                
                VORGSIZE = size(tmp);
                Ds       = VORGSIZE(2:end);
                PDs      = prod(Ds);
                APDs     = (1:PDs)';
                tmp       = reshape(tmp,[VORGSIZE(1), PDs]);
                XExt     = {cast(z3grid(:),'double'),APDs};
                F        = griddedInterpolant(XExt,tmp,'linear');
                tmp1     = cast(reshape(F({cast(zp3i(:),class(XExt{1})),APDs}),[length(zp3i), Ds]),superiorfloat(z3grid,tmp,zp3i));
                
                % Take expectation by summing over Fconditional(\tilde
                % Z_{2,t+1}, \tilde Z_{3,t+1})*prob(\tilde
                % Z_{2,t+1}, \tilde Z_{3,t+1})
                Expected                     = max(real(sum(sum(exp(tmp1).*probtilde))),10^(-320));
                return
            end
        end
                function asset = simulateStructuralIRF(asset_pI, macro_dyn, num_set)
            Nsim    = num_set.Nsim;
            Nbonds  = num_set.Nbonds;
            T       = num_set.Tirf;
            N       = num_set.N;
            sizexm  = num_set.sizexm;
            shat    = num_set.shat;
            xmgrid  = num_set.xmgrid;
            N2      = num_set.N2;
            burn    = num_set.burn;
            z       = num_set.Z;
            delta     = macro_dyn.delta;
            theta0    = macro_dyn.theta0;
            theta1    = macro_dyn.theta1;
            theta2    = macro_dyn.theta2;
            gamma     = macro_dyn.gamma;
            P         = macro_dyn.P;
            Q         = macro_dyn.Q;
            A         = macro_dyn.A;
            Ainv      = macro_dyn.Ainv;
            Sigmau    = macro_dyn.Sigmau;
            Ptilde    = macro_dyn.Ptilde;
            g         = macro_dyn.g;
            phi       = macro_dyn.phi;
            sigmac    = macro_dyn.sigmac;
            rf        = macro_dyn.rf;
            Sbar      = macro_dyn.Sbar;
            sigma_vec = macro_dyn.sigma_vec;
            
            % Copy asset prices
            asset = asset_pI;
            asset.initialShockVec = [0,0,1,0];

            % Define matrices used often
            Z14                = zeros(1,4);
            Z4T                = zeros(4,T);
            Z3T                = zeros(3,T);
            ZT1                = zeros(T,1);
            Initial            = zeros(Nsim,1);
            
            % Initialize all quantities before loop
            
            % Reshape price-consumption ratio
            G5dim              = log(reshape(asset_pI.G, N, N, N,N2,sizexm));
            G5dim_rn           = log(reshape(asset_pI.G_rn, N, N, N,N2,sizexm));
            
            % Reshape nominal bond prices
            Pnomplus5dim       = log(reshape(asset_pI.Bnom(end,:,:,:), N ,N ,N, N2, sizexm));
            Pnomminus5dim      = log(reshape(asset_pI.Bnom(end-1,:,:,:), N ,N ,N, N2, sizexm));
            Pnomplus5dim_rn    = log(reshape(asset_pI.Bnom_rn(end,:,:,:), N ,N ,N, N2, sizexm));
            Pnomminus5dim_rn   = log(reshape(asset_pI.Bnom_rn(end-1,:,:,:), N ,N ,N, N2, sizexm));
            
            % Reshape real bond prices
            Pplus5dim          = log(reshape(asset_pI.Bn(end,:,:,:), N, N, N, N2, sizexm));
            Pminus5dim         = log(reshape(asset_pI.Bn(end-1,:,:,:), N, N, N, N2, sizexm));
            Pplus5dim_rn       = log(reshape(asset_pI.Bn_rn(end,:,:,:), N, N, N, N2, sizexm));
            Pminus5dim_rn      = log(reshape(asset_pI.Bn_rn(end-1,:,:,:), N, N, N, N2, sizexm));
            
            % Reshape grid for n-dimensional interpolation
            z1shape            = reshape(z(:,1), N, N, N);
            z2shape            = reshape(z(:,2), N, N, N);
            z3shape            = reshape(z(:,3), N, N, N);
            z1grid             = reshape(z1shape(:,1,1), N, 1);
            z2grid             = reshape(z2shape(1,:,1), N, 1);
            z3grid             = reshape(z3shape(1,1,:), N, 1);
            
            % Upper and lower bounds of grids
            zgrid_             = [z1grid, z2grid, z3grid];
            zupper             = max(zgrid_)';
            zlower             = min(zgrid_)';
            supper             = max(shat);
            slower             = min(shat);
            xminusupper        = max(xmgrid);
            xminuslower        = min(xmgrid);
            
            xIrfSum = zeros(1,20);
            shatsimIrfSum = zeros (20,1); 
            cIrfSum = zeros(20,1); 
            piIrfSum = zeros(1,20);
            iIrfSum = zeros(1,20);
            PDIrfSum = zeros(1,20);
            y10nomIrfSum = zeros(1,20);
            y10realIrfSum = zeros(1,20);
            PDIrfSum_rn = zeros(1,20);
            y10nomIrfSum_rn = zeros(1,20);
            y10realIrfSum_rn = zeros(1,20);
            PDIrfSum_rp = zeros(1,20);
            y10nomIrfSum_rp = zeros(1,20);
            y10realIrfSum_rp = zeros(1,20);
            responsesSum = zeros(7,1);
            
            const2 =  [1,0,0]*( P- phi*eye(3))* Ainv;
            const3 = ([0,0,1]-[0,1,0]* P)* Ainv;
            
            %% start Nsim independent simulations
            % j counts the independent simulations
            j=1;
            while j<=Nsim
                %% Draw shocks
                % Generate T draws of u_t and compute \epsilon_t from that
                u = mvnrnd(Z14, diag((sigma_vec(1:4))), T)';
                % Run irfs in test mode without background noise
                if asset.testirf==1
                    u = Z4T;
                end
                
                % Shock period
                u(:, burn+2)   = u(:, burn+2)+(asset.initialShockVec.*sqrt(sigma_vec(1:4)))';
                
                eps            = A*Q*u;
                uast           = u(4,:);
                
                % Initialize simulated time series
                ztildesim      = Z3T;
                shatsim        = ZT1;
                Yhat           = Z3T;
                PD             = ZT1;
                PD_rn          = ZT1;
                Pnomplus       = ZT1;
                Pnomplus_rn    = ZT1;
                Pnomminus      = ZT1;
                Pnomminus_rn   = ZT1;
                Pplus          = ZT1;
                Pplus_rn       = ZT1;
                Pminus         = ZT1;
                Pminus_rn      = ZT1;
                csim           = ZT1;
                csim_ss        = ZT1; 
                piast          = ZT1;
                
                % Update state vector
                for t=3:T
                    % Dynamics for \tilde Z
                    ztildesim(:,t) = Ptilde*ztildesim(:,t-1)+eps(:,t);
                    % Dynamics for \hat Y
                    Yhat(:,t)      = P*Yhat(:,t-1)+Ainv*eps(:,t);
                    % Dynamics for surplus consumption ratio relative to steady-state
                    shatsim(t)     = ...
                        theta0*shatsim(t-1)+((1/gamma)*const3-const2)*ztildesim(:,t-1)+...
                        senshat(shatsim(t-1), Sbar)*sigmac*eps(1,t);

                    % Truncate state variables at upper and lower end of grid, so we can use standard interpolation
                    zinterp        = max(zlower, min(ztildesim(:,t), zupper));
                    sinterp        = max(slower, min(shatsim(t), supper));
                    xminusinterp   = max(xminuslower, min(const3*ztildesim(:,t-1),xminusupper));
                    csim(t)        = g + csim(t-1) + (Yhat(1,t)-phi*Yhat(1,t-1));
                    csim_ss(t) = g + csim_ss(t-1); %XXXFL
                    piast(t)       = piast(t-1)+uast(t);
                    
                    % Scaling factors for nominal bonds
                    ePplus  = exp(Nbonds*piast(t));
                    ePminus = exp((Nbonds-1)*piast(t));
                    
                    % Interpolate to obtain price-consumption ratio and real and nominal bond prices with risk premia
                    qNDInterp = eqinterpnND2(z1grid, z2grid, z3grid, shat, xmgrid, ...
                        {G5dim, Pnomplus5dim, Pnomminus5dim, Pplus5dim, Pminus5dim}, ...
                        zinterp(1), zinterp(2), zinterp(3), ...
                        sinterp, xminusinterp, 0);
                    
                    % Scale price-consumption ratio for consumption claim
                    PD(t)          = qNDInterp(1)/4;
                    
                    % 20-quarter bond prices
                    Pnomplus(t)    = qNDInterp(2)/ePplus;
                    Pplus(t)       = qNDInterp(4);
                    
                    % 19-quarter bond prices
                    Pnomminus(t)   = qNDInterp(3)/ePminus;
                    Pminus(t)      = qNDInterp(5);
                    
                    % Interpolate to obtain risk neutral price-consumption ratio and real and nominal bond prices
                    qNDInterp_rn = eqinterpnND2(z1grid, z2grid, z3grid, shat, xmgrid, ...
                        {G5dim_rn, Pnomplus5dim_rn, Pnomminus5dim_rn, Pplus5dim_rn, Pminus5dim_rn}, ...
                        zinterp(1), zinterp(2), zinterp(3), ...
                        sinterp, xminusinterp, 0);
                    
                    % Scale risk-neutral price-consumption ratio for consumption claim
                    PD_rn(t)          = qNDInterp_rn(1)/4;
                    
                    % Risk-neutral 20-quarter bond prices
                    Pnomplus_rn(t)    = qNDInterp_rn(2)/ePplus;
                    Pplus_rn(t)       = qNDInterp_rn(4);
                    
                    % Risk-neutral 19-quarter bond prices
                    Pminus_rn(t)      = qNDInterp_rn(5);
                    Pnomminus_rn(t)   = qNDInterp_rn(3)/ePminus;
                    
                end
                % 5-year nominal and real bond yields in annualized percent (in logs!)
                y10nom              = -400*log(Pnomplus)/Nbonds;
                y10real             = -400*log(Pplus)/Nbonds;
                y10nom_rn           = -400*log(Pnomplus_rn)/Nbonds;
                y10real_rn          = -400*log(Pplus_rn)/Nbonds;
                
                % Real short rate (i minus expected inflation)
                rfr                = 400*([0,0,1]*Yhat-[0,1,0]*P*Yhat+rf);
                % 1-quarter nominal yield
                rfr_nom            = 400*([0,0,1]*Yhat+ piast'+rf);
                
                x                  = 100*Yhat(1,:);
                pi                 = 400*(Yhat(2,:)+piast');
                i                  = 400*(Yhat(3,:)+piast');
                % State vector with unit root component in inflation and interest rate
                Y                  = [x;pi;i];
                
                % Levered dividends: Compute D^{delta}_{t+1}
                dDelta = PD(2:end).*exp(csim(2:end)) + exp(csim(2:end)) - (1-delta).*PD(1:end-1).*exp(csim(1:end-1)).*exp(rfr(1:end-1)./400)' ...
                    - delta.*PD(2:end).*exp(csim(2:end));
                % Take the 64-quarter moving average of D^{delta}_{t+1}
                dDelta_bar = conv(dDelta,ones(1,64),'valid')/64;
                
                % Drop burn period observations
                PD                 = PD(burn-1:end);
                PD_rn              = PD_rn(burn-1:end);
                Pnomminus          = Pnomminus(burn-1:end);
                Pnomplus           = Pnomplus(burn-1:end);
                Pminus             = Pminus(burn-1:end);
                Pplus              = Pplus(burn-1:end);
                Pnomminus_rn       = Pnomminus_rn(burn-1:end);
                Pnomplus_rn        = Pnomplus_rn(burn-1:end);
                Pminus_rn          = Pminus_rn(burn-1:end);
                Pplus_rn           = Pplus_rn(burn-1:end);
                Yhat               = Yhat(:,burn-1:end);
                Y                  = Y(:,burn-1:end);
                y10nom             = y10nom(burn-1:end);
                y10real            = y10real(burn-1:end);
                y10nom_rn          = y10nom_rn(burn-1:end);
                y10real_rn         = y10real_rn(burn-1:end);
                rfr                = rfr(burn-1:end);
                rfr_nom            = rfr_nom(burn-1:end);
                csim               = csim(burn-1:end);
                csim_ss            = csim_ss(burn-1:end); 
                dDelta             = dDelta(burn-2:end);
                dDelta_bar         = dDelta_bar((burn-65):end); 

                
                % Price-dividend ratio of levered equity at time t divides by smoothed dividends
                PDlev              = delta.*PD.*exp(csim)./dDelta_bar;
                PDlev_rn           = delta.*PD_rn.*exp(csim)./dDelta_bar;
                
                % Level return on consumption claim
                Ret              = exp((csim(2:end)-csim(1:end-1))).*(1+4*PD(2:end))./(4*PD(1:end-1));
                Ret_rn           = exp((csim(2:end)-csim(1:end-1))).*(1+4*PD_rn(2:end))./(4*PD_rn(1:end-1));
                
                % Log return on levered equity
                reteq              = log((1/delta)*Ret - ((1-delta)/delta)*exp(rfr(1:end-1)'/400));
                reteq_rn           = log((1/delta)*Ret_rn - ((1-delta)/delta)*exp(rfr(1:end-1)'/400));
                
                % In percentage units in excess of the riskfree rate
                reteq              = 100*reteq - rfr(1:end-1)'/4;
                reteq_rn           = 100*reteq_rn - rfr(1:end-1)'/4;
                
                % Log excess nominal bond returns in percent but not annualized
                retnom             = 100*log(Pnomminus(2:end))-100*log(Pnomplus(1:end-1))-rfr_nom(1:end-1)'/4;
                retnom_rn          = 100*log(Pnomminus_rn(2:end))-100*log(Pnomplus_rn(1:end-1))-rfr_nom(1:end-1)'/4;
                
                % Log excess real bond returns in percent but not annualized
                retreal            = 100*log(Pminus(2:end))-100*log(Pplus(1:end-1))-rfr(1:end-1)'/4;
                retreal_rn         = 100*log(Pminus_rn(2:end))-100*log(Pplus_rn(1:end-1))-rfr(1:end-1)'/4;
                
                % Drop observations to ensure that returns and prices have the same length
                PDlev              = PDlev(3:end);
                PDlev_rn           = PDlev_rn(3:end);
                reteq              = reteq(2:end-1);
                reteq_rn           = reteq_rn(2:end-1);
                retnom             = retnom(2:end-1);
                retnom_rn          = retnom_rn(2:end-1);
                retreal            = retreal(2:end-1);
                retreal_rn         = retreal_rn(2:end-1);
                Yhat               = Yhat(:,3:end);
                Y                  = Y(:,3:end);
                y10nom             = y10nom(3:end);
                y10real            = y10real(3:end);
                y10nom_rn          = y10nom_rn(3:end);
                y10real_rn         = y10real_rn(3:end);
                rfr                = rfr(3:end);
                rfr_nom            = rfr_nom(3:end);
                csim               = csim(3:end);
                csim_ss            = csim_ss(3:end); 
                
                if min(PDlev)<0
                    j=j+1;
                    continue
                end
                % Increase loop counter
                j = j+1;
                
                xIrf       = Y(1,1:20) - Y(1,1);
                cIrf       = (csim(1:20) - csim_ss(1:20)) - (csim(1) - csim_ss(1)); 
                piIrf      = Y(2,1:20) - Y(2,1);
                iIrf       = Y(3,1:20) - Y(3,1);
                PDIrf      = (cumsum(reteq(1:20)-reteq(1)));
                y10nomIrf  = y10nom(1:20) - y10nom(1);
                y10realIrf = y10real(1:20) - y10real(1);
                shatsimIrf = shatsim(burn+1:burn+20);
                
                responses = [reteq(2); % Return on stocks on impact
                    reteq_rn(2);
                    (y10nom(2) - y10real(2)) - (y10nom(1) - y10real(1)); % Change in breakeven
                    (y10nom_rn(2) - y10real_rn(2)) - (y10nom_rn(1) - y10real_rn(1)); % Change in breakeven
                    rfr_nom(2) - rfr_nom(1); % change in short rate
                    Y(1,2) - Y(1,1); % Change in output on impact
                    Y(1,8) - Y(1,1)]; % Change in output in 8 quarters
                
                % Risk-neutral vs. Risk premium
                PDIrf_rn       = (cumsum(reteq_rn(1:20)-reteq_rn(1)));
                y10nomIrf_rn   = y10nom_rn(1:20) - y10nom_rn(1);
                y10realIrf_rn  = y10real_rn(1:20) - y10real_rn(1);
                PDIrf_rp       = PDIrf - PDIrf_rn;
                y10nomIrf_rp   = y10nomIrf - y10nomIrf_rn;
                y10realIrf_rp  = y10realIrf - y10realIrf_rn;
                
                xIrfSum = xIrfSum + xIrf;
                shatsimIrfSum = shatsimIrfSum + shatsimIrf; 
                cIrfSum = cIrfSum + cIrf; 
                piIrfSum = piIrfSum + piIrf;
                iIrfSum = iIrfSum + iIrf;
                PDIrfSum = PDIrfSum + PDIrf';
                y10nomIrfSum = y10nomIrfSum +y10nomIrf';
                y10realIrfSum = y10realIrfSum + y10realIrf';
                
                PDIrfSum_rn = PDIrfSum_rn + PDIrf_rn';
                y10nomIrfSum_rn = y10nomIrfSum_rn +y10nomIrf_rn';
                y10realIrfSum_rn = y10realIrfSum_rn + y10realIrf_rn';
                
                PDIrfSum_rp = PDIrfSum_rp + PDIrf_rp';
                y10nomIrfSum_rp = y10nomIrfSum_rp +y10nomIrf_rp';
                y10realIrfSum_rp = y10realIrfSum_rp + y10realIrf_rp';
                responsesSum = responsesSum + responses;
            end
            xIrf = xIrfSum/Nsim;
            cIrf = cIrfSum/Nsim; 
            shatsimIrf = shatsimIrfSum/Nsim; 
            piIrf = piIrfSum/Nsim;
            iIrf = iIrfSum/Nsim;
            PDIrf = PDIrfSum/Nsim;
            y10nomIrf = y10nomIrfSum/Nsim;
            y10realIrf = y10realIrfSum/Nsim;
            
            PDIrf_rn = PDIrfSum_rn/Nsim;
            y10nomIrf_rn = y10nomIrfSum_rn/Nsim;
            y10realIrf_rn = y10realIrfSum_rn/Nsim;
            
            PDIrf_rp = PDIrfSum_rp/Nsim;
            y10nomIrf_rp = y10nomIrfSum_rp/Nsim;
            y10realIrf_rp = y10realIrfSum_rp/Nsim;
            responses = responsesSum/Nsim;
            
            responsesObj.eqRet = responses(1);
            responsesObj.eqRet_rn = responses(2);
            responsesObj.breakevenChange = responses(3);
            responsesObj.breakevenChange_rn = responses(4);
            responsesObj.ffrNomChange = responses(5);
            responsesObj.outputChange = responses(6);
            responsesObj.outputIn8QrtChange = responses(7);
            
            Irftemp.x=xIrf;
            Irftemp.c = cIrf; 
            Irftemp.pi=piIrf;
            Irftemp.i=iIrf;
            Irftemp.PD=PDIrf;
            Irftemp.shatsim=shatsimIrf;
            Irftemp.y10nom=y10nomIrf;
            Irftemp.y10real=y10realIrf;
            Irftemp.PD_rn=PDIrf_rn;
            Irftemp.y10nom_rn=y10nomIrf_rn;
            Irftemp.y10real_rn=y10realIrf_rn;
            Irftemp.PD_rp=PDIrf_rp;
            Irftemp.y10nom_rp=y10nomIrf_rp;
            Irftemp.y10real_rp=y10realIrf_rp;
            Irftemp.responses = responsesObj;
            
            asset.Irftemp    =   Irftemp;
            return
        end
    end
end
