clear variables
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure A5: Scatter Bond Betas by Shocks and Monetary Policy %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load asset moment results, name can be Supply (for Panel A), Demand (for Panel B) and MP (for Panel C). All the results consider zeta=0.0 
name='Supply';
input_file = ['./results_sensitivity/loop_results/results_gridsearch_',name,'.mat'];
load(input_file)

% Filter to keep only reasonable values for gamma_x (between 0.6 and 0.8).
a=and(Table1All(12,:)>=0.6, Table1All(12,:)<=0.8);
Table1All=Table1All(:,a);
Table3All=Table3All(:,a);

% Plot bond betas against MP inflation weight
h=figure;
scatter(Table1All(10,:), Table3All(10,:), 'o') % Scatter plot for nominal bond beta
hold on
scatter(Table1All(10,:), Table3All(11,:),'x') % Scatter plot for real bond beta
yline(0)
axis([1.1,2,-0.3,0.4])
xlabel('MP Inflation Weight \gamma^\pi','fontsize',14)
ylabel('Bond-Stock Betas','fontsize',14)
if strcmp(name, 'MP')
    legend('Nominal Bond Beta','Real Bond Beta','Location','SouthWest')
    legend boxoff
end
hold off
saveas(h,['./figures/',name,'_counterfactuals_gamma_pi_Pflueger.png'],'png'); % Save the figure as a PNG file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Table 3: Model-Implied Post-Pandemic Monetary Policy and Shocks %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Best match for the bond-stock betas for each period. You must consider
% name='Supply' to replicate Table 3 results

% Condition on particular parameter values, such as rho_i = 0.8 and gamma^x=0.5
a=and(Table1All(11,:)==0.5,Table1All(12,:)==0.8);
Table1All=Table1All(:,a);
Table3All=Table3All(:,a);

% Target values for each row (Nominal and Real betas)
target_values = [-0.01, 0.05; 0.49, 0.39]; % Specify target values for nominal and real betas
data_close = []; % Initialize a variable to store the closest match data
diff=zeros(1,max(size(target_values))); % Initialize a variable to store the closest match data

for ix=1:size(target_values,1)
    % Calculate the absolute difference from the target values for each row
    abs_diff_row1 = abs(Table3All(10,:) - target_values(ix,1)).^2;
    abs_diff_row2 = abs(Table3All(11,:) - target_values(ix,2)).^2;
    
    % Combine the differences 
    combined_diff = abs_diff_row1 + abs_diff_row2;
    
    % Find the index of the minimum combined difference
    [~, columnIndex] = min(combined_diff);
    
    % Display the column index of the closest value
    disp(['Column index of the closest values to the targets of nominal and real betas (', num2str(target_values(ix,1)), ' and ', ...
        num2str(target_values(ix,2)), ') for period ', num2str(ix),' : ', num2str(columnIndex)]);
    
    % Prepare the data for the table
    data_close = [Table1All(10:12,columnIndex); Table1All(14:16,columnIndex); Table3All(10:11,columnIndex)];    
    
    if ix==1
        data_close_all = data_close; % Initialize the result table on the first iteration
    else
        data_close_all = [data_close_all data_close]; % Append results for subsequent iterations
    end
    diff(ix)=combined_diff(columnIndex); % Store the minimum difference for the current target values
end

% Create a table to display the results
results=vertcat(target_values(:,1)', target_values(:,2)',[data_close_all(end-1:end,1), data_close_all(end-1:end,2)], [data_close_all(1:end-2,1), data_close_all(1:end-2,2)], diff);
results=round(results,2);
resultTable = table(results);
resultTable.Properties.RowNames = {'Data Nominal', 'Data Real', 'Nominal Bond-Stock Beta', 'Real Bond-Stock Beta','\gamma^\pi', '\gamma^x', '\rho^i', '\sigma_x', '\sigma_\pi', '\sigma_ST',  append('Distance',name)};

% Display the table
resultTable
