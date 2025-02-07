%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation of asset pricing moments for Table 3 and Figure A5 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop for gridsearch changing rho^i, gamma^pi and gamma^x with a MP shock
clear variables

% Change directory to the specified path
cd('C:\Users\luisyepezsa\OneDrive - The University of Chicago\Documents\GitHub\ProgrammingPackage_public\Pflueger1980s_replication_files');

% Load baseline calibration
load Calibrationgridsearch_demand2.mat

% Set seed for replication
rng(0);

% Loop for gridsearch
for ii=0:0.1:0.7 % Set grid values for gamma_x
    for jj=0:0.1:0.9 % Set grid values for gamma_pi
        for kk=0:0.1:0.4 % Set grid values for rho_i
            % Monetary policy rule
            macro.gamma_x          =   (0.4+ii)/4;        % MP Response to output (annualized percent)
            macro.gamma_pi         =   1.1+jj;            % MP Response to inflation
            macro.rho_i            =   0.5+kk;            % MP Persistence
            
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
            
            % Adaptive inflation expectations parameter chosen to match Campbell-Shiller
            % regressions
            macro.zeta=0.0;
            
            % Variance of quarterly MP shock in natural units
            macro.sigma_vec(3)     =   (1.19/400)^2;   
            macro.sigma_vec        = [10^(-9) 10^(-9) macro.sigma_vec(3) 10^(-9)];  
           
            % Solve the model
            [Table1, Table3]    =   generateTable3(macro, num);

            if ii==0 && jj==0 && kk==0
                Table3All = Table3;
                Table1All = Table1;
            else
                Table3All = [Table3All Table3];
                Table1All = [Table1All Table1];
            end  

            
            % Save the result (change the path and format as needed)
            save('./results_sensitivity/loop_results/results_gridsearch_MP.mat',"Table1All","Table3All");
        end
    end
end

%% Loop for gridsearch changing rho^i, gamma^pi and gamma^x with a Supply shock
clear variables

% Change directory to the specified path
cd('C:\Users\luisyepezsa\OneDrive - The University of Chicago\Documents\GitHub\ProgrammingPackage_public\Pflueger1980s_replication_files');

% Load baseline calibration
load Calibrationgridsearch_demand2.mat

% Set seed for replication
rng(0);

% Loop for gridsearch
for ii=0:0.1:0.7 % Set grid values for gamma_x
    for jj=0:0.1:0.9 % Set grid values for gamma_pi
        for kk=0:0.1:0.4 % Set grid values for rho_i
            % Monetary policy rule
            macro.gamma_x          =   (0.4+ii)/4;        % MP Response to output (annualized percent)
            macro.gamma_pi         =   1.1+jj;            % MP Response to inflation
            macro.rho_i            =   0.5+kk;            % MP Persistence
            
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
            
            % Adaptive inflation expectations parameter chosen to match Campbell-Shiller
            % regressions
            macro.zeta             =   0.0;
            
            % Variance of quarterly shocks in natural units 
            macro.sigma_vec(1)     =   (0.10/100)^2; 
            macro.sigma_vec(2)     =   (0.58/400)^2; 
            macro.sigma_vec(3)     =   (0.55/400)^2;   
            macro.sigma_vec        = [macro.sigma_vec(1) macro.sigma_vec(2)  macro.sigma_vec(3) 10^(-9)];   
           
            % Solve the model
            [Table1, Table3]    =   generateTable3(macro, num);

            if ii==0 && jj==0 && kk==0
                Table3All = Table3;
                Table1All = Table1;
            else
                Table3All = [Table3All Table3];
                Table1All = [Table1All Table1];
            end  

            
            % Save the result (change the path and format as needed)
            save('./results_sensitivity/loop_results/results_gridsearch_Supply.mat',"Table1All","Table3All");
        end
    end
end

%% Loop for gridsearch changing rho^i, gamma^pi and gamma^x with a Demand shock
clear variables

% Change directory to the specified path
cd('C:\Users\luisyepezsa\OneDrive - The University of Chicago\Documents\GitHub\ProgrammingPackage_public\Pflueger1980s_replication_files');

% Load baseline calibration
load Calibrationgridsearch_demand2.mat

% Set seed for replication
rng(0);

% Loop for gridsearch
for ii=0:0.1:0.7 % Set grid values for gamma_x
    for jj=0:0.1:0.9 % Set grid values for gamma_pi
        for kk=0:0.1:0.4 % Set grid values for rho_i
            % Monetary policy rule
            macro.gamma_x          =   (0.4+ii)/4;        % MP Response to output (annualized percent)
            macro.gamma_pi         =   1.1+jj;            % MP Response to inflation
            macro.rho_i            =   0.5+kk;            % MP Persistence
            
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
            
            % Adaptive inflation expectations parameter chosen to match Campbell-Shiller
            % regressions
            macro.zeta             =   0.0;
            
            % Variance of quarterly shocks in natural units 
            macro.sigma_vec(1)     =   (0.59/100)^2; 
            macro.sigma_vec(2)     =   (0.07/400)^2; 
            macro.sigma_vec(3)     =   (0.07/400)^2;   
            macro.sigma_vec        = [macro.sigma_vec(1) macro.sigma_vec(2)  macro.sigma_vec(3) 10^(-9)];   
           
            % Solve the model
            [Table1, Table3]    =   generateTable3(macro, num);

            if ii==0 && jj==0 && kk==0
                Table3All = Table3;
                Table1All = Table1;
            else
                Table3All = [Table3All Table3];
                Table1All = [Table1All Table1];
            end  

            
            % Save the result (change the path and format as needed)
            save('./results_sensitivity/loop_results/results_gridsearch_Demand.mat',"Table1All","Table3All");
        end
    end
end
