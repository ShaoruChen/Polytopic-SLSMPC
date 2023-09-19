% plot the coverage of each robust MPC method with varying eps_A parameters
% for horizon = 3 and 10. 

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
fig = tiledlayout(1,2, 'TileSpacing','Compact');
set(gcf, 'Position', [100, 100, 1000, 750]);


 nexttile
 hold on 
for ii = 1:4
    plot(eps_A_list, best_coverage_mat(:,ii), 'o-.','MarkerSize', 14, 'LineWidth', 3.0 );
end

for ii = 5:7
    plot(eps_A_list, best_coverage_mat(:,ii), 's-', 'MarkerSize', 14, 'LineWidth', 3.0 );
end


grid on
set(gca, 'XTick',0.05:0.05:0.45);
ax = gca; % get handle to the current axes
ax.FontSize = 22;
text(0.75, 0.85, 'hor.=3', 'FontSize', 24, 'FontWeight', 'bold','Units', 'normalized')

ylim([0 1]);
xlabel('$\epsilon_A$', 'Interpreter', 'Latex', 'FontSize', 28);
ylabel('coverage', 'Interpreter', 'Latex', 'FontSize', 28);


%% load data of hor = 10

data = load('data/example_2d_diags_list_eps_A.mat');
diags_list_eps_A = data.diags_list;
eps_A_list = 0.05:0.05:0.45;

N_eps = length(eps_A_list);
method_cell = {'Tube_MPC', 'Tube_MPC_Homothetic', 'Tube_MPC_Nominal', 'Tube_MPC_Flexible',  ...
               'SLS_MPC', 'Lumped_Dist_MPC', 'Constr_Tightening_MPC'};
legend_cell = {'Tube-A', 'Tube-B', 'Tube-C','Tube-D','SLS-MPC', 'Lumped-Disturbance', 'Offline-Tightening' }
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

ax = nexttile;
hold on 
for ii = 1:4
    h(ii) = plot(eps_A_list, best_coverage_mat(:,ii), 'o-.','MarkerSize', 14, 'LineWidth', 3.0 );
    % set(h(ii), {'DisplayName'}, legend_cell{ii});
    h(ii).DisplayName = legend_cell{ii};
end

for ii = 5:7
    h(ii) = plot(eps_A_list, best_coverage_mat(:,ii), 's-', 'MarkerSize', 14, 'LineWidth', 3.0 );
    % set(h(ii), {'DisplayName'}, legend_cell{ii});
        h(ii).DisplayName = legend_cell{ii};

end


grid on
set(gca, 'XTick',0.05:0.05:0.45);
ax = gca; % get handle to the current axes
ax.FontSize = 22;

text(0.75, 0.85, 'hor.=10', 'FontSize', 24, 'FontWeight', 'bold','Units', 'normalized')

ylim([0 1]);
xlabel('$\epsilon_A$', 'Interpreter', 'Latex', 'FontSize', 28);
ylabel('coverage', 'Interpreter', 'Latex', 'FontSize', 28);

%% add legend

lg = legend(ax, h(:), 'Orientation', 'Horizontal', 'NumColumns', 4);
lg.Layout.Tile = 'North';

exportgraphics(fig,'eps_A_hor.png','BackgroundColor','none');

