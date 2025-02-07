%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation of asset pricing moments for Figures 8 and A6 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop for gridsearch with the second period calibration (demand shock) and changing rho^i, gamma_pi and sigma_pi - Figure 8 results
clear variables

% Change directory to the specified path
cd('C:\Users\luisyepezsa\OneDrive - The University of Chicago\Documents\GitHub\InflationExpectationsFormation');

% Load demand shock calibration
load Calibrationgridsearch_demand2.mat

% Set seed for replication
rng(0);

% Set grid values for rho^i
grid_rho_i =[0.54 0.80];

% Loop for gridsearch
for ii=grid_rho_i
    for jj=0:0.1:0.9 % Set grid values for gamma_pi
        for kk=0:1:12 % Set grid values for sigma_pi
            % Monetary policy rule
            macro.gamma_x          =   1/4;            % MP Response to output (annualized percent)
            macro.gamma_pi         =   1.1+jj;         % MP Response to inflation
            macro.rho_i            =   ii;             % MP Persistence
                    
            % Preference parameters
            macro.theta0           =   0.9658;         % Peristence surplus consumption
            macro.theta1           =   -0.84;          % Backward looking habit
            macro.phi              =   0.9900;         % Consumption-output gap                                
            macro.gamma            =   2;              % Utility curvature
            macro.g                =   0.004725;       % Consumption growth      
            macro.rf               =   0.00235;        % Risk-free rate
                    
            % Phillips curve
            macro.kappa            =   0.0062/4;       % Slope of the Phillips curve
                   
            % Leverage parameter
            macro.delta            =   0.66; 
                    
            %adaptive inflation expectations parameter chosen to match Campbell-Shiller
            %regressions
            macro.zeta             =   0.0;
                    
            % Variance of quarterly shocks in natural units. 
            macro.sigma_vec(1)     =   (0.59/100)^2; 
            macro.sigma_vec(2)     =   ((0.07*(kk+1))/400)^2;  
            macro.sigma_vec(3)     =   (0.07/400)^2;   
            macro.sigma_vec        = [macro.sigma_vec(1) macro.sigma_vec(2)  macro.sigma_vec(3) 10^(-9)];   
                   
            % Solve the model
            [Table1, Table3]    =   generateTable3(macro, num);
        
            if ii==0.54 && jj==0 && kk==0
                Table3All = Table3;
                Table1All = Table1;
            else
                Table3All = [Table3All Table3];
                Table1All = [Table1All Table1];
            end  
                    
            % Save the result (change the path and format as needed)
            save('./results_sensitivity/loop_results/results_gridsearch_second_period_rho_i_gamma_pi_sigma_pi.mat',"Table1All","Table3All");
        end
    end
end
%% Loop for gridsearch with the second period calibration (demand shock) and changing rho^i, gamma^x and sigma_pi - Figure A6 results
clear variables

% Change directory to the specified path
cd('C:\Users\luisyepezsa\OneDrive - The University of Chicago\Documents\GitHub\InflationExpectationsFormation');

% Load demand shock calibration
load Calibrationgridsearch_demand2.mat

% Set seed for replication
rng(0);

% Set grid values for rho^i
grid_rho_i =[0.54 0.80];

% Loop for gridsearch
for ii=grid_rho_i
    for jj=0:0.1:0.7 % Set grid values for gamma_pi
        for kk=0:1:12 % Set grid values for sigma_pi
            % Monetary policy rule
            macro.gamma_x          =   (0.4+jj)/4;     % MP Response to output (annualized percent)
            macro.gamma_pi         =   1.1;            % MP Response to inflation
            macro.rho_i            =   ii;             % MP Persistence
                    
            % Preference parameters
            macro.theta0           =   0.9658;         % Peristence surplus consumption
            macro.theta1           =   -0.84;          % Backward looking habit
            macro.phi              =   0.9900;         % Consumption-output gap                                
            macro.gamma            =   2;              % Utility curvature
            macro.g                =   0.004725;       % Consumption growth      
            macro.rf               =   0.00235;        % Risk-free rate
                    
            % Phillips curve
            macro.kappa            =   0.0062/4;       % Slope of the Phillips curve
                   
            % Leverage parameter
            macro.delta            =   0.66; 
                    
            %adaptive inflation expectations parameter chosen to match Campbell-Shiller
            %regressions
            macro.zeta             =   0.0;
                    
            % Variance of quarterly shocks in natural units. 
            macro.sigma_vec(1)     =   (0.59/100)^2; 
            macro.sigma_vec(2)     =   ((0.07*(kk+1))/400)^2;  
            macro.sigma_vec(3)     =   (0.07/400)^2;   
            macro.sigma_vec        = [macro.sigma_vec(1) macro.sigma_vec(2)  macro.sigma_vec(3) 10^(-9)];   
                   
            % Solve the model
            [Table1, Table3]    =   generateTable3(macro, num);
        
            if ii==0.54 && jj==0 && kk==0
                Table3All = Table3;
                Table1All = Table1;
            else
                Table3All = [Table3All Table3];
                Table1All = [Table1All Table1];
            end  
                    
            % Save the result (change the path and format as needed)
            save('./results_sensitivity/loop_results/results_gridsearch_second_period_rho_i_gamma_x_sigma_pi.mat',"Table1All","Table3All");
        end
    end
end