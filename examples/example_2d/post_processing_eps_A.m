clear
close all


%% eps_A comparison
data = load('data/example_2d_diags_list_eps_A_hor_3.mat');
diags_list_eps_A = data.diags_list;
eps_A_list = 0.05:0.05:0.45;

N_eps = length(eps_A_list);
method_cell = {'Tube_MPC', 'Tube_MPC_Homothetic', 'Tube_MPC_Nominal', 'Tube_MPC_Flexible',  ...
               'SLS_MPC', 'Lumped_Dist_MPC', 'Constr_Tightening_MPC'};
num_method = 7;

best_coverage_mat = zeros(N_eps, num_method);
solver_time_mat = zeros(N_eps, num_method);
for ii = 1:N_eps
    diags_record = diags_list_eps_A{ii};
    N = length(diags_record);
    
    % find the best coverage for each tube-based MPC method
    best_coverage = zeros(1, length(method_cell));
    for jj = 1:N
        diags = diags_record{jj};
        method = diags.method;
        ind = find(strcmp(method_cell, method));
        
        coverage = diags.feasible_rate;
        if coverage >= best_coverage(ind)
            best_coverage(ind) = coverage;
            solver_time_mat(ii, ind) = diags.avg_runtime_feasible;
        end
    end
    best_coverage_mat(ii, :) = best_coverage;
end
figure(1);
set(gcf, 'Position', [100, 100, 1000, 750]);


subplot(121)
hold on 
for ii = 1:4
    plot(eps_A_list, best_coverage_mat(:,ii), 'o-.','MarkerSize', 14, 'LineWidth', 3.0 );
end

for ii = 5:7
    plot(eps_A_list, best_coverage_mat(:,ii), 's-', 'MarkerSize', 14, 'LineWidth', 3.0 );
end

ylim([0 1]);
xlabel('$\epsilon_A$', 'Interpreter', 'Latex', 'FontSize', 32);
ylabel('coverage', 'Interpreter', 'Latex', 'FontSize', 32);
grid on
set(gca, 'XTick',0.05:0.05:0.45);
ax = gca; % get handle to the current axes
ax.FontSize = 22;
text(0.75, 0.85, 'hor.=3', 'FontSize', 24, 'FontWeight', 'bold','Units', 'normalized')


%% load data of hor = 10

data = load('data/example_2d_diags_list_eps_A.mat');
diags_list_eps_A = data.diags_list;
eps_A_list = 0.05:0.05:0.45;

N_eps = length(eps_A_list);
method_cell = {'Tube_MPC', 'Tube_MPC_Homothetic', 'Tube_MPC_Nominal', 'Tube_MPC_Flexible',  ...
               'SLS_MPC', 'Lumped_Dist_MPC', 'Constr_Tightening_MPC'};
num_method = 7;

best_coverage_mat = zeros(N_eps, num_method);
solver_time_mat = zeros(N_eps, num_method);
for ii = 1:N_eps
    diags_record = diags_list_eps_A{ii};
    N = length(diags_record);
    
    % find the best coverage for each tube-based MPC method
    best_coverage = zeros(1, length(method_cell));
    for jj = 1:N
        diags = diags_record{jj};
        method = diags.method;
        ind = find(strcmp(method_cell, method));
        
        coverage = diags.feasible_rate;
        if coverage >= best_coverage(ind)
            best_coverage(ind) = coverage;
            solver_time_mat(ii, ind) = diags.avg_runtime_feasible;
        end
    end
    best_coverage_mat(ii, :) = best_coverage;
end

subplot(122)
hold on 
for ii = 1:4
    h(ii) = plot(eps_A_list, best_coverage_mat(:,ii), 'o-.','MarkerSize', 14, 'LineWidth', 3.0 );
end

for ii = 5:7
    h(ii) = plot(eps_A_list, best_coverage_mat(:,ii), 's-', 'MarkerSize', 14, 'LineWidth', 3.0 );
end

ylim([0 1]);
xlabel('$\epsilon_A$', 'Interpreter', 'Latex', 'FontSize', 32);
ylabel('coverage', 'Interpreter', 'Latex', 'FontSize', 32);
grid on
set(gca, 'XTick',0.05:0.05:0.45);
ax = gca; % get handle to the current axes
ax.FontSize = 22;

text(0.75, 0.85, 'hor.=10', 'FontSize', 24, 'FontWeight', 'bold','Units', 'normalized')

%% add legend


% Create legends
lgd1 = legend([h(1), h(2), h(3), h(4)], 'Tube-A', 'Tube-B', 'Tube-C', 'Tube-D', 'NumColumns', 4, 'Interpreter', 'latex', 'FontSize', 24);
lgd1.Box = 'off';
axes_copy = copyobj(gca, gcf);

lgd2 = legend(axes_copy, [h(5), h(6), h(7)], 'SLS-MPC', 'Lumped-Disturbance', 'Offline-Tightening', 'NumColumns', 3, 'Interpreter', 'latex', 'FontSize', 24);
lgd2.Box = 'off';

% Get the positions
pos1 = get(lgd1, 'Position');
pos2 = get(lgd2, 'Position');

% Adjust the position of the second legend
set(lgd2, 'Position', [pos1(1) + pos1(3) - pos2(3), pos1(2), pos2(3), pos2(4)]);

% Hide the copied axes
set(axes_copy, 'Visible', 'off');


% % Create legends
% legend1 = legend([h(1), h(2), h(3), h(4)], 'First', 'Second', 'Third', 'Fourth', 'NumColumns', 4, 'Location', 'southoutside');
% legend1.Box = 'off';
% 
% % Make a copy of the axes, which will contain the second legend
% axes_copy = copyobj(gca, gcf);
% 
% % Create the second legend
% legend2 = legend(axes_copy, [h(5), h(6), h(7)], 'Fifth', 'Sixth', 'Seventh', 'NumColumns', 3, 'Location', 'southoutside');
% legend2.Box = 'off';
% 
% % Get the positions
% pos1 = get(legend1, 'Position');
% pos2 = get(legend2, 'Position');
% 
% % Adjust the position of the second legend
% set(legend2, 'Position', [pos1(1) + pos1(3) - pos2(3), pos1(2), pos2(3), pos2(4)]);
% 
% % Hide the copied axes
% set(axes_copy, 'Visible', 'off');


% Lgnd = legend('Tube-A', 'Tube-B', 'Tube-C', 'Tube-D', ...
%                'SLS-MPC', 'Lumped-Disturbance', 'Offline-Tightening', 'NumColumns', 4, ...
%                'Interpreter', 'latex', 'FontSize', 14, 'Location', 'bestoutside');

% Lgnd = legend('show');
% Lgnd.Position(1) = 0.01;
% Lgnd.Position(2) = -0.01;

%% set figure size
% Adjust the figure properties
% set(gcf, 'PaperPositionMode', 'auto')
% 
% % Save the figure as a PNG file with no border
% print('myplot.png', '-dpng', '-r300')

