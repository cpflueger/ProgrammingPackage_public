clear variables
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 8: 3D Nominal Bond-Stock Beta Scatter Plots for rho^i, gamma^pi and sigma_pi values %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the data file
input_file = './results_sensitivity/loop_results/results_gridsearch_second_period_rho_i_gamma_pi_sigma_pi.mat';
load(input_file)

% Values to loop over
rho_i_values = [0.54 0.80];

% Loop through each value of rho^i
for idx = 1:length(rho_i_values)
    rho_i = rho_i_values(idx);

    % Filter the data for the current rho^i
    filter_idx = Table1All(12,:) == rho_i;
    
    % Define your data for axes based on the filter
    x = Table1All(15, filter_idx); % sigma_pi values (x-axis)
    y = Table1All(10, filter_idx); % gamma^pi values (y-axis)
    z_nominal = Table3All(10, filter_idx); % Nominal Bond Beta (z-axis)
    z_real = Table3All(11, filter_idx); % Real Bond Beta (z-axis)

    % Create grid for target planes
    [xq, yq] = meshgrid(linspace(min(x), max(x), 50), linspace(min(y), max(y), 50));

    % 3D Scatter Plot for Nominal Bond Beta
    h1 = figure;
    scatter3(x, y, z_nominal, 50, 'filled', 'MarkerEdgeColor', 'b', 'HandleVisibility', 'off'); % Scatter plot with blue color for nominal bond betas
    hold on;

    % Add a target plane at z = 0.24 with contrasting color
    target_nominal = 0.24 * ones(size(xq));
    surf(xq, yq, target_nominal, 'FaceAlpha', 0.7, 'FaceColor', [1, 0.8, 0.6], 'EdgeColor', 'none', 'DisplayName',''); % Orange plane for target with legend

    % Customize plot appearance
    ax = gca;
    ax.FontSize = 12; 
    xlabel('Volatility Supply Shock \sigma_\pi', 'fontsize', 12)
    ylabel({'MP Inflation', 'Weight \gamma^\pi'}, 'fontsize', 12)
    zlabel('Nominal Bond-Stock Beta', 'fontsize', 12)
    grid on
    view(3) % 3D view
    axis([min(x), max(x), min(y), max(y), -0.35, 0.3])
    
    % Save the nominal scatter plot
    saveas(h1, ['./figures/scatter_gamma_pi_Pflueger_', num2str(rho_i*100)], 'png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure A6: 3D Nominal Bond-Stock Beta Scatter Plots for rho^i, gamma^x and sigma_pi values %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the data file
input_file = './results_sensitivity/loop_results/results_gridsearch_second_period_rho_i_gamma_x_sigma_pi.mat';
load(input_file)

% Values to loop over
rho_i_values = [0.54 0.80];

% Loop through each value of rho^i
for idx = 1:length(rho_i_values)
    rho_i = rho_i_values(idx);

    % Filter the data for the current rho^i
    filter_idx = Table1All(12,:) == rho_i;
    
    % Define your data for axes based on the filter
    x = Table1All(15, filter_idx); % sigma_pi values (x-axis)
    y = Table1All(11, filter_idx); % gamma^x values (y-axis)
    z_nominal = Table3All(10, filter_idx); % Nominal Bond Beta (z-axis)
    z_real = Table3All(11, filter_idx); % Real Bond Beta (z-axis)

    % Create grid for target planes
    [xq, yq] = meshgrid(linspace(min(x), max(x), 50), linspace(min(y), max(y), 50));

    % 3D Scatter Plot for Nominal Bond Beta
    h1 = figure;
    scatter3(x, y, z_nominal, 50, 'filled', 'MarkerEdgeColor', 'b', 'HandleVisibility', 'off'); % Scatter plot with blue color for nominal bond betas
    hold on;

    % Add a target plane at z = 0.24 with contrasting color
    target_nominal = 0.24 * ones(size(xq));
    surf(xq, yq, target_nominal, 'FaceAlpha', 0.7, 'FaceColor', [1, 0.8, 0.6], 'EdgeColor', 'none', 'DisplayName', ''); % Orange plane for target with legend

    % Customize plot appearance
    ax = gca;
    ax.FontSize = 12; 
    xlabel('Volatility Supply Shock \sigma_\pi', 'FontSize', 12) 
    ylabel({'MP Output', 'Weight \gamma^x'}, 'FontSize', 12) 
    zlabel('Nominal Bond-Stock Beta', 'fontsize', 12)
    grid on
    view(3) % 3D view
    axis([min(x), max(x), min(y), max(y), -0.35, 0.3])

    % Save the nominal scatter plot
    saveas(h1, ['./figures/scatter_gamma_x_Pflueger_', num2str(rho_i*100)], 'png');
end

