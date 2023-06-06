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
                                   %    R-squared 1-year excess return onto
                                   %    pd
        nominalBonds,              % Structure of moments related to nominal bonds:
                                   %    Term premium, bond return volatility, Sharpe ratio, mean log yield spread, 
                                   %    volatility log yield spread, AR(1) coefficient log yield spread, 
                                   %    coefficient 1-year bond excess returns onto lagged log yield spread,
                                   %    R^2 1-year bond excess returns onto lagged log yield spread
        realBonds,                 % Structure of moments related to real bonds:
                                   %    Term premium, volatility of real bond returns, Sharpe ratio, mean log yield spread,                                
                                   %    volatility log yield spread, real bond-stock beta, correlation real bond returns with stock returns
        crossAsset,                % Structure of cross-asset moments:
                                   %    correlation nominal bond returns
                                   %    with stock returns
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
              
        AP_responses,              %   coefficients of bond and stock returns components (cf, real rate, risk premium) onto impulse responses
        CovDecomposition,          %   covariances of bond and stock return components (cf, real rate, risk premium)
        CorrDecomposition          %   correlations of bond and stock return components (cf, real rate, risk premium)
        additional_moments,        %   fraction simulation periods with s_t>s^max, std approximation error for 1-period nominal rate
        simulated,                 %   simulated series that can reproduce some important moments 
        simulated_rn,              %   simulated series that can reproduce some important moments when we are in the risk neutral case
    end
%% Methods section
    methods
        %% Initialize output asset_p class and save implied parameters to solve for asset prices
        function asset = initial_ap(assetInput,macro_dyn)
        % Initialize output asset_p class
        asset = asset_p(assetInput.G,assetInput.Bn,assetInput.Bnom,0);
        % Set up for risk premia version (that we use then when evaluate
        % the asset price recursion)
        asset.risk_neutral_run = 0; 
        % Save implied parameters
        asset.ImpliedParams    = macro_dyn.ImpliedParams;
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
                % Set a preliminary value for log-linearization constant.
                % It does not enter into any numbers shown in the paper and
                % Simply serves to avoid an error message
                assetRN.rho = 0.99;
                % Simulate risk neutral prices and returns
                disp('Simulate risk neutral moments')
                tic
                assetRN = assetRN.SimulateMoments(num_set,macro_dyn);
                toc
                % Risk-neutral price-dividend ratio
                PDRN = assetRN.stocks.meanPDlev;
                % Log-linearization constant to decompose risk-neutral
                % Returns into cash flow news and discount rate news
                asset.rho = 1/(1+1/(4*PDRN));
                % Save risk-neutral asset prices in aset
                asset.G_rn = assetRN.G;
                asset.Bn_rn = assetRN.Bn;
                asset.Bnom_rn = assetRN.Bnom;
                asset.simulated_rn = assetRN.simulated;
            else
                asset = assetInput;
                asset.ImpliedParams    = macro_dyn.ImpliedParams;
                asset.G       = zeros(num_set.N1,num_set.N2,num_set.sizexm);
                asset.Bn      = zeros(num_set.Nbonds-1,num_set.N1,num_set.N2,num_set.sizexm);
                asset.Bnom    = zeros(num_set.Nbonds-1,num_set.N1,num_set.N2,num_set.sizexm);
                asset.G_rn    = zeros(num_set.N1,num_set.N2,num_set.sizexm);
                asset.Bn_rn   = zeros(num_set.Nbonds-1,num_set.N1,num_set.N2,num_set.sizexm);
                asset.Bnom_rn = zeros(num_set.Nbonds-1,num_set.N1,num_set.N2,num_set.sizexm);
                asset.rho     = 0.9956; 
            end
        end
        %% Constructor method
        function asset = asset_p(G , Bn, Bnom, ImpliedParams)         
        % this method initializes the class given some input arguments
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
%           asset.CorrDecomposition = correlation matrix of stock and bond
%           cash flow news, real rate news, and risk premium return components
%           asset.CovDecomposition = covariance matrix corresponding to asset.CorrDecomposition
%           asset.additional moments = fraction of surplus consumption
%           ratio simulations above smax, std of approximation error for
%           one-period nominal interest rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function asset  = SimulateMoments(asset_pI, num_set, macro_dyn)
            Nsim    = num_set.Nsim;
            T       = num_set.T;
            N       = num_set.N;
            sizexm  = num_set.sizexm;
            shat    = num_set.shat;
            xmgrid  = num_set.xmgrid;
            N2      = num_set.N2;
            burn    = num_set.burn;
            z       = num_set.Z;
            gamma     = macro_dyn.gamma;
            delta     = macro_dyn.delta;
            theta0    = macro_dyn.theta0;
            theta1    = macro_dyn.theta1; 
            theta2    = macro_dyn.theta2;
            P         = macro_dyn.P;
            Q         = macro_dyn.Q;
            Ainv      = macro_dyn.Ainv;
            Sigmau    = macro_dyn.Sigmau;
            Ptilde    = macro_dyn.Ptilde;
            g         = macro_dyn.g;
            phi       = macro_dyn.phi;
            sigmac    = macro_dyn.sigmac;
            vast      = macro_dyn.vast;
            rf        = macro_dyn.rf;
            sigmaperp = macro_dyn.sigmaperp;
            Sbar      = macro_dyn.Sbar;
            smax      = macro_dyn.smax; 
            QM        =[1,0,0]*macro_dyn.Q;

            % Copy asset prices
            asset = asset_p(asset_pI.G,asset_pI.Bn,asset_pI.Bnom,macro_dyn.ImpliedParams);          
            asset.G_rn    = asset_pI.G_rn;
            asset.Bn_rn   = asset_pI.Bn_rn;
            asset.Bnom_rn = asset_pI.Bnom_rn;
            asset.rho     = asset_pI.rho;
            asset.simulated_rn     = asset_pI.simulated_rn;
            
            % Define matrices used often 
            Z13                = zeros(1,3);
            Z3T                = zeros(3,T);
            ZT1                = zeros(T,1);
            Initial            = zeros(Nsim,1);

            %Initialize all quantities before loop
            rfrNomStd          = Initial;
            rho_nom_rfr        = Initial;
            std_dp             = Initial;
            mean_pdlev         = Initial;
            rho_dp             = Initial;
            re1                = Initial;
            re1_r2             = Initial;
            ys1                = Initial;
            ys1_r2             = Initial;
            TermSlope          = Initial;
            TermSlopeStd       = Initial;
            TermSlopeReal      = Initial;
            TermSlopeRealStd   = Initial;
            EqPremium          = Initial;
            BondPremium        = Initial;
            RealBondPremium    = Initial;
            fracmax            = Initial;
            AR_slope5          = Initial;
            consGrowthVol      = Initial;
            consGrowthAR1      = Initial;
            xVol               = Initial;
            xAR1               = Initial;
            piChangeVol        = Initial;
            piChangeAR1        = Initial;
            correlations       = repmat(Initial, 1,2);
            Stdeq              = Initial;
            Stdnom             = Initial;
            Stdreal            = Initial;
            beta_nom           = Initial';
            beta_real          = Initial';
            stdApproxErrI      = Initial;
            CovDecomposition   = zeros(Nsim, 6, 6);
            CorrDecomposition  = zeros(Nsim, 6, 6);
            

            % Reshape price-consumption ratio
            G5dim              = log(reshape(asset_pI.G, N, N, N,N2,sizexm));
            G5dim_rn           = log(reshape(asset_pI.G_rn, N, N, N,N2,sizexm));
            
            % Reshape nominal bond prices 
            Pnomplus5dim       = log(reshape(asset_pI.Bnom(end,:,:,:), N ,N ,N, N2, sizexm));
            Pnomminus5dim      = log(reshape(asset_pI.Bnom(end-1,:,:,:), N ,N ,N, N2, sizexm)); 
            Pnomplus5dim_rn       = log(reshape(asset_pI.Bnom_rn(end,:,:,:), N ,N ,N, N2, sizexm));
            Pnomminus5dim_rn      = log(reshape(asset_pI.Bnom_rn(end-1,:,:,:), N ,N ,N, N2, sizexm)); 
            
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
            
            %% start Nsim independent simulations
            % j counts the independent simulations
           j=1;
           while j<=Nsim 
                %% Draw shocks
                % Generate T draws of \epsilon_t
                eps            = mvnrnd(Z13, diag([1,1,0]), T)';

                %Independent shock to unit root component of inflation
                epsperp        = mvnrnd(ZT1, 1, T)';
                %shock to unit root component of inflation=component driven
                %by \epsilon_t + independent shock
                uast           = vast'*eps+sigmaperp*epsperp;
                
                %initialize simulated time series
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
                piast          = ZT1;
                approxErrI     = ZT1;

                %update state vector
                for t=3:T
                    % dynamics for \tilde Z
                    ztildesim(:,t) = Ptilde*ztildesim(:,t-1)+eps(:,t);
                    % dynamics for \hat Y
                    Yhat(:,t)      = P*Yhat(:,t-1)+Ainv*eps(:,t);
                    % dynamics for surplus consumption ratio relative to
                    % steady-state
                    shatsim(t)     = ...
                        theta0*shatsim(t-1)+theta1*Yhat(1,t-1)+theta2*Yhat(1,t-2)+...
                        senshat(shatsim(t-1), Sbar)*sigmac*eps(1,t);

                    %truncate state variables at upper and lower end of
                    %grid, so we can use standard interpolation
                    zinterp        = max(zlower, min(ztildesim(:,t), zupper));
                    sinterp        = max(slower, min(shatsim(t), supper));
                    xminusinterp   = max(xminuslower, min(Yhat(1,t-1),xminusupper));
                    csim(t)        = g+csim(t-1)+(Yhat(1,t)-phi*Yhat(1,t-1));
                    piast(t)       = piast(t-1)+uast(t);
                     
                    % scaling factors for nominal bonds
                    ePplus  = exp(20*piast(t));
                    ePminus = exp(19*piast(t));
                    
                    %interpolate to obtain price-consumption ratio and real
                    %and nominal bond prices with risk premia
                    qNDInterp = eqinterpnND2(z1grid, z2grid, z3grid, shat, xmgrid, ...
                                            {G5dim, Pnomplus5dim, Pnomminus5dim, Pplus5dim, Pminus5dim}, ...
                                            zinterp(1), zinterp(2), zinterp(3), ...
                                            sinterp, xminusinterp, 0);
                    
                    %scale price-consumption ratio for consumption claim                    
                    PD(t)          = qNDInterp(1)/4;
                    
                    %20-quarter bond prices
                    Pnomplus(t)    = qNDInterp(2)/ePplus;
                    Pplus(t)       = qNDInterp(4);

                    %19-quarter bond prices
                    Pnomminus(t)   = qNDInterp(3)/ePminus;
                    Pminus(t)      = qNDInterp(5);
                    
                    %interpolate to obtain risk neutral price-consumption ratio and real
                    %and nominal bond prices 
                    qNDInterp_rn = eqinterpnND2(z1grid, z2grid, z3grid, shat, xmgrid, ...
                                            {G5dim_rn, Pnomplus5dim_rn, Pnomminus5dim_rn, Pplus5dim_rn, Pminus5dim_rn}, ...
                                            zinterp(1), zinterp(2), zinterp(3), ...
                                            sinterp, xminusinterp, 0);
                    
                    %scale risk-neutral price-consumption ratio for
                    %consumption claim 
                    PD_rn(t)          = qNDInterp_rn(1)/4;
                    
                    %risk-neutral 20-quarter bond prices
                    Pnomplus_rn(t)    = qNDInterp_rn(2)/ePplus;
                    Pplus_rn(t)       = qNDInterp_rn(4);
                    
                    %risk-neutral 19-quarter bond prices
                    Pminus_rn(t)      = qNDInterp_rn(5);
                    Pnomminus_rn(t)   = qNDInterp_rn(3)/ePminus;
                    
                    % approximation error for 1-period nominal rate
                    approxErrI(t) = .5*([0,1,0]*Q(:,2:4) + [0,0,1])*Sigmau(2:4,2:4)*([0,1,0]*Q(:,2:4) + [0,0,1])'+ gamma*(senshat(shatsim(t), Sbar) + 1)*QM(2:4)*Sigmau(2:4,2:4)*([0,1,0]*Q(:,2:4) + [0,0,1])';     
                end
                % 5-year nominal and real bond yields in annualized percent (in logs!)
                y5nom              = -400*log(Pnomplus)/20;
                y5real             = -400*log(Pplus)/20;
                
                % Real short rate (i minus expected inflation)
                rfr                = 400*([0,0,1]*Yhat-[0,1,0]*P*Yhat+rf);
                
                % Levered dividends: Compute D^{delta}_{t+1} 
                dDelta = PD(2:end).*exp(csim(2:end)) + exp(csim(2:end)) - (1-delta).*PD(1:end-1).*exp(csim(1:end-1)).*exp(rfr(1:end-1)./400)' - delta.*PD(2:end).*exp(csim(2:end));
                % Take the 64-quarter moving average of D^{delta}_{t+1}
                dDelta_bar = conv(dDelta,ones(1,64),'valid')/64;
                                
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
                y5nom              = y5nom(burn-1+2:end); 
                y5real             = y5real(burn-1+2:end); 
                rfr                = rfr(burn-1+2:end);
                dDelta_bar         = dDelta_bar((burn-65)+2:end); 
                eps                = eps(:,burn-1+3:end);
                approxErrI         = approxErrI(burn-1+2:end); 
                burn1              = burn+2; 
                
                % Price-dividend ratio of levered equity at time t divides
                % by smoothed dividends
                PDlev              = delta.*PD.*exp(csim)./dDelta_bar;	
                
                %1-quarter nominal yield
                rfr_nom            = 400*([0,0,1]*Yhat+piast+rf);
                
                % Level return on consumption claim
                Ret              = exp((csim(2:end)-csim(1:end-1))).*(1+4*PD(2:end))./(4*PD(1:end-1));
                Ret_rn           = exp((csim(2:end)-csim(1:end-1))).*(1+4*PD_rn(2:end))./(4*PD_rn(1:end-1));

                % log return on levered equity
                reteq              = log((1/delta)*Ret - ((1-delta)/delta)*exp(rfr(1:end-1)'/400));
                reteq_rn           = log((1/delta)*Ret_rn - ((1-delta)/delta)*exp(rfr(1:end-1)'/400));

                % in percentage units in excess of the riskfree rate
                reteq              = 100*reteq - rfr(1:end-1)'/4;
                reteq_rn           = 100*reteq_rn - rfr(1:end-1)'/4;
                
                % log excess nominal bond returns in percent but not
                % annualized
                retnom             = 100*log(Pnomminus(2:end))-100*log(Pnomplus(1:end-1))-rfr_nom(1:end-1)'/4;
                retnom_rn          = 100*log(Pnomminus_rn(2:end))-100*log(Pnomplus_rn(1:end-1))-rfr_nom(1:end-1)'/4;
                
                % log excess real bond returns in percent but not
                % annualized
                retreal            = 100*log(Pminus(2:end))-100*log(Pplus(1:end-1))-rfr(1:end-1)'/4;
                retreal_rn         = 100*log(Pminus_rn(2:end))-100*log(Pplus_rn(1:end-1))-rfr(1:end-1)'/4;
                               
                % Drop observations to ensure that returns and prices have
                % the same length
                
                %nominal and real log yield spreads
                spreadNom          	= y5nom'-rfr_nom;
                spreadReal          = y5real'-rfr;
                
                % 1-year log equity excess returns in natural units
                ret1yr             = conv(reteq,ones(1,4),'valid')/100;
                                
                %levered price-dividend ratio
                pdlev=log(PDlev);
                               
                % compute cash-flow news and real rate news
                % vectors to compute equity real rate news analytically according to
                % Campbell and Ammer 
                rho                = asset_pI.rho;              
                rhoIPinv           = inv(eye(3)-rho*P);
                Gammaeq_rr            = -rho*([0,0,1]-[0,1,0]*P)*rhoIPinv*Ainv;
                   
                % cash-flow news of stock and bond returns
                reteq_cf           = reteq_rn-(100*Gammaeq_rr*eps)';
                retnom_cf      = retnom_rn-retreal_rn;
                               
                % risk premium excess returns of stock and bond returns
                reteq_rp=reteq-reteq_rn;
                retnom_rp=retnom-retnom_rn;  
                
                %% compute covariance matrix and correlation matrix of stock
                %and bond return decomposition
                %covariances in percent
                CovDecomposition(j,:,:) = cov([retnom_cf, retnom_rn-retnom_cf, retnom_rp, reteq_cf, reteq_rn-reteq_cf, reteq_rp]);
                CorrDecomposition(j,:,:)= corrcoef([retnom_cf, retnom_rn-retnom_cf, retnom_rp, reteq_cf, reteq_rn-reteq_cf, reteq_rp]);
                
                %% Stock moments
                
                % Equity risk premium in annualized units
                EqPremium(j)       = 4*(mean(reteq) + .5*std(reteq)^2/100);
                %std. log equity excess returns
                Stdeq(j)           = std(reteq)*2;
                
                % exp(mean(log pd))
                mean_pdlev(j)     = exp(mean(pdlev));
                % Standard deviation of log dp
                std_dp(j)          = std(pdlev);               
                % Autocorrelation of dp
                dp_corr            = corrcoef(pdlev(2:end), pdlev(1:end-1));
                rho_dp(j)          = dp_corr(1,2);
                
                % Predictability with pd
                % 1 quarter, 1 year and 5 year regressions of returrn on price-dividend ratios
                [re1_coef,~,~,~,R2_re1]           = regress(ret1yr, [ones(size(pdlev(1:end-4))), pdlev(1:end-4)]);                               
                re1(j)             = re1_coef(2);
                re1_r2(j)          = R2_re1(1);

                %% Nominal bond moments
               
                %nominal term premium
                BondPremium(j)     = 4*(mean(retnom) + .5*std(retnom)^2/100);
                %std of bond returns
                Stdnom(j)          = std(retnom)*2;               
                
                %mean and std of log yield spread
                TermSlope(j)       = mean(spreadNom);
                TermSlopeStd(j)    = std(spreadNom);
                
                % Autocorrelation of log yield spread 
                AR_slope5temp      = corrcoef(spreadNom(2:end), spreadNom(1:end-1));
                AR_slope5(j)       = AR_slope5temp(1,2);              
                
                % Regress returns onto lagged log yield spread
                ret1yrNom             = retnom(4:end)+retnom(3:end-1)+retnom(2:end-2)+retnom(1:end-3);
                ret1yrNom             = ret1yrNom/100;               
                
                % multiply returns by 100 to match units in
                % empirical exercise     
                [ys1_coef,~,~,~,R2_ys1]           = regress(100*ret1yrNom, [ones(size(spreadNom(1:end-4)')), spreadNom(1:end-4)']);
                ys1(j)             = ys1_coef(2);
                ys1_r2(j)          = R2_ys1(1);               
                
                %% Cross-asset 
                
                %bond-stock return correlations
                bondstock_corr_temp     = corrcoef(retnom, reteq);
                tipsstock_corr_temp     = corrcoef(retreal, reteq);
                correlations(j,:)       = [bondstock_corr_temp(1,2), tipsstock_corr_temp(1,2)];
               
                %nominal bond beta
                beta_temp           = regress(retnom, [ones(T-burn1+1,1), reteq]);
                beta_nom(j)         = beta_temp(2);
                                    
                
                %% Real bonds
                
                %Term premium
                RealBondPremium(j) = 4*(mean(retreal) + .5*std(retreal)^2/100);
                %std returns
                Stdreal(j)         = std(retreal)*2;
                
                %Mean and std. log yield spread
                TermSlopeReal(j)   = mean(spreadReal);
                TermSlopeRealStd(j)= std(spreadReal);
                % Real bond beta
                beta_temp          = regress(retreal, [ones(T-burn1+1,1), reteq]);
                beta_real(j)       = beta_temp(2);          
                
                
                %% macro dynamics
                
                %std and AR(1) of changes in nominal 1-quarter yield
                rfrNomStd(j)       = std(rfr_nom(2:end)-rfr_nom(1:end-1));
                rfr_nom_corr       = corrcoef(rfr_nom(3:end)-rfr_nom(2:end-1),rfr_nom(2:end-1)-rfr_nom(1:end-2));
                rho_nom_rfr(j)     = rfr_nom_corr(1,2);
                
                %inflation changes
                piChanges            = 4*(100*Yhat(2,1:end-1)+100*piast(1:end-1) - (100*Yhat(2,2:end)+100*piast(2:end)));
                piChangeVol(j)       = std(piChanges);
                piChangeAR1temp      = corrcoef(piChanges(1:end-1), piChanges(2:end));
                piChangeAR1(j)       = piChangeAR1temp(1,2);
                
                %log consumption growth
                consGrowthVol(j)     = 2*std(100*(csim(2:end)-csim(1:end-1)));
                consGrowthAR1temp    = corrcoef(csim(2:end-1)-csim(1:end-2),csim(3:end)-csim(2:end-1));
                consGrowthAR1(j)     = consGrowthAR1temp(1,2);                      
                
                %output gap
                xVol(j)              = std(100*Yhat(1,:));
                xAR1temp             = corrcoef(Yhat(1,2:end), Yhat(1,1:end-1));
                xAR1(j)              = xAR1temp(1,2);
          
                %fraction s_t>s^max
                fracmax(j)         = sum(shatsim+log(Sbar)>smax)/T;
                %std of approximation error for 1-period nominal rate in
                %annualized percent
                stdApproxErrI(j)     = std(approxErrI)*400;
                
                %skip simulation run if PDlev turns negative
                if min(PDlev)<0
                    j=j+1;
                    continue
                end
                %increase loop counter
                j=j+1;               
            end
            % End of loop average across Nsim simulations
            
            %% save output 
            asset.CovDecomposition = reshape(mean(CovDecomposition),6,1);
            asset.CorrDecomposition = reshape(mean(CorrDecomposition),6,1);
            
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
            nominalBonds.coeffRegRetOnYS1y         = ys1;
            nominalBonds.R2RegRetOnYS1y            = ys1_r2;

            asset.nominalBonds                     = nominalBonds;
     
            %% Real Bonds
            realBonds.termPremium          = RealBondPremium;
            realBonds.vol                  = mean(Stdreal);
            realBonds.sharpeRatio          = RealBondPremium/mean(Stdreal);
            realBonds.meanLogYieldSpread   = TermSlopeReal;
            realBonds.volLogYieldSpread    = TermSlopeRealStd;
            realBonds.betaRealStock        = beta_real;
            realBonds.corrRealStock        = correlations(2);
            
            asset.realBonds                = realBonds;
            
            %% Cross-Asset moments
            crossAsset.corrNomStock     = correlations(1);
            crossAsset.betaNom          = beta_nom;
            
            asset.crossAsset            = crossAsset;
            
            %% Simulated series

            simulated.pdlev     = pdlev;
            simulated.Yhat      = Yhat;
            simulated.retnom    = retnom;
            simulated.reteq     = reteq;
            simulated.y5nom     = y5nom;
            simulated.rfr_nom   = rfr_nom;  
            simulated.retreal   = retreal;
            simulated.piast     = piast;
            
            asset.simulated     = simulated;
            %% Macro Dynamics moments
            macroDynamics.iChangeVol     = rfrNomStd; 
            macroDynamics.iChangeAR1     = rho_nom_rfr;  
            
            macroDynamics.piChangeVol    = mean(piChangeVol);
            macroDynamics.piChangeAR1    = mean(piChangeAR1);
            
            macroDynamics.consGrowthVol  = mean(consGrowthVol);
            macroDynamics.consGrowthAR1  = mean(consGrowthAR1);
            
            macroDynamics.xVol           = mean(xVol);
            macroDynamics.xAR1           = mean(xAR1);
            
            asset.macroDynamics          = macroDynamics;
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
            xGL3 = 0*num_set.xGL3;
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

            % Initialize matrices for log asset prices bnew=log(Bn),
            % bnnew=log(Bnom), fnew=log(Fn)=log price-consumption ratio of zero-coupon
            % consumption claim
            bnew  = zeros(N1,N2, sizexm);
            bnnew = zeros( N1, N2, sizexm);
            fnew  = zeros( N1, N2, sizexm);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Evaluate two-period bond prices and one-period
            %consumption-claim using analytic expressions
            
            %Implement analytic expressions for one-period zero-coupon
            %consumption claim 
            
            %asset prices with risk premia
            if risk_neutral_run == 0
                % Define temporary constants outside the loop to calculate them only once
                const1 =  g;
                const2 =  [1,0,0]*( P- phi*eye(3))* Ainv;
                const3 = ([0,0,1]-[0,1,0]* P)* Ainv;
                %i1 loops over grid points for \tilde Z
                for i1=1: N1
                    % i2 loops over grid points for \hat s_t
                    for i2=1: N2
                        % i3 loops over grid points for x_{t-1}
                        for i3=1: sizexm
                            fnew(i1,i2,i3)= const1 + const2*Z(i1,:)' - rf - const3*Z(i1,:)' - 0.5*gamma*(1- theta0)*(1-2* shat(i2)) + 0.5*(gamma* lambdas(i2)+ (gamma- 1))^2*sigmac^2;
                        end
                    end
                end
                
            %risk-neutral asset prices   
            else
                % Define temporary constants for risk-neutral asset prices outside the loop to calculate them only once
                const1 =  g;
                const2 =  [1,0,0]*( P- phi*eye(3))* Ainv;
                const3 = ([0,0,1]-[0,1,0]* P)* Ainv;
                %i1 loops over grid points for \tilde Z
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
            
            %define temporary constants outside the loop to calculate them only once
            const1 = -([0,0,1]-[0,1,0]* P)*(eye(3)+ P)* Ainv;
            const2 = -[0,0,1]*(eye(3)+ P)* Ainv;
            vr2     = ([0,0,1]-[0,1,0]* P)* Q;
            vbn2    = [0,1,1]*Q+2*[0,0,0,1];
            
            %asset priceswith risk premia
            if risk_neutral_run == 0
                %i1 loops over grid points for \tilde Z
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
                
            %risk-neutral asset prices    
            else
                %i1 loops over grid points for \tilde Z
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
            %prepare inputs to speed up the evaluation of the
            %numerical expectation outside the loop that evaluates the
            %recursion
            
            % distribution of \tilde Z_{t+1} conditional on \tilde Z_t       
            % initialize distributions for three dimensions of \tilde
            % Z_{t+1}
            zp               = zeros(N1, GLpoints);
            zp2              = zeros(N1, GLpoints2);
            zp3              = zeros(N1, GLpoints3);
            % i1 loops over grid for \tilde Z_t
            for i1=1: N1
                %E \tilde Z_{t+1} conditional on Z_t 
                ztmp         =  Ptilde* Z(i1,:)';
                %Distribution of first dimension of \tilde Z_{t+1} conditional on \tilde Z_t
                zp(i1,:)     = [1,0,0]*ztmp+ xGL;
                %Distribution of second dimension of \tilde Z_{t+1} conditional on \tilde Z_t
                zp2(i1,:)    = [0,1,0]*ztmp+ xGL2;
                %Distribution of first dimension of \tilde Z_{t+1} conditional on \tilde Z_t. 
                % Notice that zp3 is constant across its third dimension, because the shock \epsilon_{3,t+1} has zero variance by construction.                
                zp3(i1,:)    = [0,0,1]*ztmp+ 0*xGL3;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % create map between Z and Ztilde. 
            
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
            % prep grid and probabilities so we can take expectations over
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

            % split up vec*\epsilon_{t+1}=vec*(1)\epsilon_{1,t+1}
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
            %inflation shocks driven by \epsilon_{2,t+1} and
            %\epsilon_{3,t+1}. vpi_mesh varies with \epsilon_{2,t+1} along rows and with \epsilon_{3,t+1} along columns  
            vpi2         = vpi(2)* xGL2';
            vpi3         = vpi(3)* xGL3';
            [vpia2, vpia3] = meshgrid(vpi2, vpi3);
            vpi_mesh       = vpia2+vpia3;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % prep interpolation along grid for x^- to improve speed
            % within the main recursion
            
            %for each value in X(i1) weightlow(i1) is such that
            %X(i1)=weightlow(i1)*xmgrid(1)+(1-weightlow(i1)*xmgrid(2)
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
            %initialize matrices used to extract asset prices along
            %sub-dimensions of the grid within the loop to help with speed and 
            %memory
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
            % evaluate asset pricing recursion
            % n loops over the maturity of zero coupon bonds/consumption claims
            % The maturity of the zero-coupon consumption claim equals n and the maturity of the 
            % zero coupon bonds equals n+1
            
            % A note on notation within the recursion: Variables ending on _old
            % refer to (n-1)-period zero coupon consumption claims and n-period zero
            % coupon bonds. Variables ending on _plus refer to distributions
            % at time t+1.
            for n=2: NN
                % log of n-1 period price-dividend ratio, real bond price and nominal bond price
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

                        % i3 loops over lagged output gap x_{t-1}
                        for i3=1: sizexm
                            % distribution of time t+1 log surplus
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
                                % find points in Z corresponding to
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
                                    
                                    %interpolate f_{n-1,t+1}(\tilde
                                    %Z_{1,t+1},s_{t+1},x_t) over \tilde
                                    %Z_{1,t+1}, i.e. first dimension of
                                    %\tilde Z_{t+1}
                                    FS1         = (1-T1).*diag(fmesh_old(XIPOS2, XIPOS1))+T1.*diag(fmesh_old(XIPOS2, XIPOS1+1));
                                    FS2         = (1-T1).*diag(fmesh_old(XIPOS2+1, XIPOS1))+T1.*diag(fmesh_old(XIPOS2+1, XIPOS1+1));
                                                                        
                                    %interpolate f_{n-1,t+1}(\tilde
                                    %Z_{1,t+1},s_{t+1},x_t) over s_{t+1}
                                    fs_plus          = (1-T2').*FS1+T2'.*FS2;
                                    
                                    %this is the recursion equation for
                                    %consumption claims as detailed in
                                    %Appendix section D.1.3
                                    if risk_neutral_run == 0
                                        fs_plus = fs_plus' + const2 + const1 * Z(i1,:)' - rf - const3 * Z(i1,:)' + const4 * (1 - 2* shat(i2)) - (gamma * (1+ lambdas(i2))- 1) * ve1;
                                    else 
                                        fs_plus = fs_plus' + const2 + const1 * Z(i1,:)' - rf - const3 * Z(i1,:)' + ve1;
                                    end
                                Fconditional(i4)      = real(exp(fs_plus)* prob1);
                                
                                %for n<= maximum bond maturity do the same
                                %calucation to obtain Bconditional and
                                %Bnconditional
                                if n < Nbonds
                                    % Interpolate log-linearly to evaluate
                                    % b_{n-1,t+1}(\tilde
                                    % Z_{t+1},s_{t+1},x_t) and b^\$_{n-1,t+1}(\tilde
                                    % Z_{t+1},s_{t+1},x_t) over x_t
                                    bmesh_old             = reshape(bn_old(a>0,:,mapxlow(i1)),  N,  N2)'*weightlow(i1)+reshape(bn_old(a>0,:,mapxhigh(i1)),  N,  N2)'*(1-weightlow(i1));
                                    Bnmesh_old            = reshape(bnom_old(a>0,:,mapxlow(i1)), N,  N2)'*weightlow(i1)+reshape(bnom_old(a>0,:,mapxhigh(i1)),  N,  N2)'*(1-weightlow(i1));
                                    
                                    % Interpolate b_{n-1,t+1}(\tilde
                                    %Z_{1,t+1},s_{t+1},x_t) over \tilde
                                    %Z_{1,t+1} and s_{t+1}
                                        FS1         = (1-T1).*diag(bmesh_old(XIPOS2, XIPOS1))+T1.*diag(bmesh_old(XIPOS2, XIPOS1+1));
                                        FS2         = (1-T1).*diag(bmesh_old(XIPOS2+1, XIPOS1))+T1.*diag(bmesh_old(XIPOS2+1, XIPOS1+1));
                                        bs_plus          = (1-T2').*FS1+T2'.*FS2;

                                    % Interpolate b^\$_{n-1,t+1}(\tilde
                                    %Z_{1,t+1},s_{t+1},x_t) over \tilde
                                    %Z_{1,t+1} and s_{t+1}
                                        FS1         = (1-T1).*diag(Bnmesh_old(XIPOS2, XIPOS1))+T1.*diag(Bnmesh_old(XIPOS2, XIPOS1+1));
                                        FS2         = (1-T1).*diag(Bnmesh_old(XIPOS2+1, XIPOS1))+T1.*diag(Bnmesh_old(XIPOS2+1, XIPOS1+1));
                                        bns_plus         = (1-T2').*FS1+T2'.*FS2;
                                        
                                        
                                    %these are the recursion equations for
                                    %real and nominal bonds as detailed in
                                    %Appendix section D.1.3    
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
                            
                            %distributions of \tilde Z_{2,t+1} and \tilde
                            %Z_{3,t+1} conditional on \tilde Z_t
                            zp2i                      = real(zp2(i1,:)');
                            zp3i                      = real(zp3(i1,:)');
                            
                            %take expectation of Fconditional over \tilde
                            %Z_{2,t+1} and \tilde Z_{3,t+1}
                            sheet(i2,i3) = expectedinterpolated(Fconditional, zp2i, zp3i, probtilde,  N, X2, X3);
                            
                            if n< Nbonds
                            %take expectation over Bconditional over \tilde
                            %Z_{2,t+1} and \tilde Z_{3,t+1}
                                sheetB(i2,i3) = expectedinterpolated(Bconditional, zp2i, zp3i, probtilde,  N, X2, X3);
                                
                            %Before taking expectations over Bnconditional
                            %multiply probabilities by the shock to bond
                            %cash flows due to inflation that is perfectly
                            %correlated with \epsilon_{2,t+1} and
                            %\epsilon_{3,t+1}
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

                %bound Fconditional away from zero to avoid error messages
                Fconditional                 = max(Fconditional, 10^(-320));
                %take log - we interpolate log-linearly
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
                %[XExt,YExt]  = ndgrid(zp2i(:),APDs);
                %finterp2       =  interpn(z2grid(:),APDs,fcond,XExt, YExt);
                %finterp2       =  cast(reshape(interpn(z2grid(:),APDs,fcond,XExt, YExt),[length(zp2i), Ds]),superiorfloat(z2grid,fcond,zp2i));
                F        = griddedInterpolant(XExt,fcond,'linear');
                finterp2       = cast(reshape(F({cast(zp2i(:),class(XExt{1})),APDs}),[length(zp2i), Ds]),superiorfloat(z2grid,fcond,zp2i));
                tmp                          = finterp2';

                VORGSIZE = size(tmp);
                Ds       = VORGSIZE(2:end);
                PDs      = prod(Ds);
                APDs     = (1:PDs)';
                tmp       = reshape(tmp,[VORGSIZE(1), PDs]);
                XExt     = {cast(z3grid(:),'double'),APDs};
                %[XExt,YExt]  = ndgrid(zp3i(:),APDs);
                %tmp1       =  interpn(z3grid(:),APDs,tmp,XExt, YExt);
                %tmp1       =  cast(reshape(interpn(z3grid(:),APDs,tmp,XExt, YExt),[length(zp3i), Ds]),superiorfloat(z3grid,tmp,zp3i));
                F        = griddedInterpolant(XExt,tmp,'linear');
                tmp1     = cast(reshape(F({cast(zp3i(:),class(XExt{1})),APDs}),[length(zp3i), Ds]),superiorfloat(z3grid,tmp,zp3i));

                %take expectation by summing over Fconditional(\tilde
                %Z_{2,t+1}, \tilde Z_{3,t+1})*prob(\tilde
                %Z_{2,t+1}, \tilde Z_{3,t+1})
                Expected                     = max(real(sum(sum(exp(tmp1).*probtilde))),10^(-320));
            return 
            end
        end
                     
    end
end