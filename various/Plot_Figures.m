%% Script to plot figures, mesh diagnostics and water budget

% Select appropriate workspace from the Analysis folder and run script to 
% produce plots of pressure head and water content at the specified 
% parameters/time. 
%
% For example, choosing workspace 'CSG 1km_5yearSteadyState - 10years.mat'
% and running the script displays the pressure head and water content
% plots at a simulation time of 10 years using a 5 year 
% Steadystate with the CSG located at 1km and switched on.
%
% When plotting water budget, select the WaterBudget workspace and toggle
% WaterBudget = true on the script.
%
% When plotting mesh diagnostics, select the Gradient/Uniform Mesh
% Diagnostics workspace and toggle Mesh_Diagnostics = true.  

% Toggle on Water Budget and Mesh Diagnostics. 
WaterBudget = false;         % Turn on when plotting Water Budget 
Mesh_Diagnostics = false;    % Turn on when plotting mesh diagnostics

% Error outputted when both are turned on
if WaterBudget == true && Mesh_Diagnostics == true
    error('Please only turn on WaterBudget or Mesh Diagnostics')
end 

% Pressure head plot
figure
surf(cfg.X, cfg.Z, vec2mat(h,r)', 'EdgeColor', 'none', 'FaceColor', 'interp')
view(2)
title('Pressure Head', 'Interpreter', 'LateX', 'FontSize', 20)
xlabel('Width (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
ylabel('Height (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
colorbar
caxis([-15 0])

% Water Content plot
figure
surf(cfg.X, cfg.Z, vec2mat(psiI,r)', 'EdgeColor', 'none', 'FaceColor', 'interp')
view(2)
title('Water Content', 'Interpreter', 'LateX', 'FontSize', 20)
xlabel('Width (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
ylabel('Height (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
colorbar
caxis([0 0.4686])

% Plots time of each iteration
if Mesh_Diagnostics == true
    % Plot the loop iteration times
    figure
    plot(1:1:iter-1, loopTimings(~isnan(loopTimings)))
    title('looptimings per iteration')
    xlabel('Iterations')
    ylabel('Loop Timings')
    sum(loopTimings(~isnan(loopTimings)));
    
    % Plot the change in step size
    figure
    title('Change in Step Size and k', 'Interpreter', 'LateX', 'FontSize', 20)
        xlabel('Iterations', 'Interpreter', 'LateX', 'FontSize', 16)
        % left axis
        yyaxis left
        plot(dt_vec, '-', 'LineWidth', 1.5)
        ylabel('Step Size (dt)', 'Interpreter', 'LateX', 'FontSize', 16)
        yyaxis right
        plot(1:length(k_vec), k_vec, '--', 'LineWidth', 1.5)
        ylim([0, 15])
        ylabel('Convergence of k', 'Interpreter', 'LateX', 'FontSize', 16)
end

% Plot Water Budget
if WaterBudget == true
    % Plot the total water flow into/out of the aquifer
    figure
    plot(t_vec, tot_vec, '--', 'LineWidth', 1.5)
    title('Water Budget', ...
            'Interpreter', 'LateX', 'FontSize', 20)
    xlabel('Time (days)', 'Interpreter', 'LateX', 'FontSize', 16)
    ylabel('Water Content', 'Interpreter', 'LateX', 'FontSize', 16)
    legend('Total water flux', ...
            'Location', 'SouthWest')
end
