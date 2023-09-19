
%% post processing

data = load('data/MPC_feasibility_random_4.mat', 'diags_list');
diags_list_4 = data.diags_list;
results_4 = diags_list_4(find(~cellfun(@isempty, diags_list_4)));

data = load('data/MPC_feasibility_random_3.mat', 'diags_list');
diags_list_3 = data.diags_list;
results_3 = diags_list_3(find(~cellfun(@isempty, diags_list_3)));

data = load('data/MPC_feasibility_random_2.mat', 'diags_list');
diags_list_2 = data.diags_list;
results_2 = diags_list_2(find(~cellfun(@isempty, diags_list_2)));

data = load('data/MPC_feasibility_random_1.mat', 'diags_list');
diags_list_1 = data.diags_list;
results_1 = diags_list_1(find(~cellfun(@isempty, diags_list_1)));

results = [results_4 results_3 results_2 results_1];
% save('feasibility_random_results.mat', 'result');


N_examples = size(results, 2);

method_cell = {'Tube_MPC', 'Tube_MPC_Homothetic', 'Tube_MPC_Nominal', 'Tube_MPC_Flexible',  ...
               'SLS_MPC', 'Lumped_Dist_MPC', 'Constr_Tightening_MPC'};
num_method = 7;

best_coverage_mat = zeros(N_examples, num_method);
solver_time_mat = zeros(N_examples, num_method);

for ii = 1:N_examples
    diags_record = results{ii};
    N = length(diags_record);

    contain_empty_diags = 0;
    for jj = 1:N
        if isempty(diags_record{jj})
            contain_empty_diags = 1;
        end
    end

    if contain_empty_diags == 1
        continue
    end

    % find the best coverage for each tube-based MPC method
    best_coverage = zeros(1, length(method_cell));
    for jj = 1:N
        diags = diags_record{jj};
        if isempty(diags)
            keyboard;
        end

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

[slsmpc_sorted, ind] = sort(best_coverage_mat(:,5));

figure;
hold on 
for ii = 1:4
    coverage = best_coverage_mat(:,ii);
    plot(coverage(ind), 's-.', 'LineWidth', 1.0 , 'LineStyle', 'none' );
end

for ii = 5:7
    coverage = best_coverage_mat(:,ii);
    if ii == 5
        plot(coverage(ind), 's-', 'LineWidth', 1.0);
    else
        plot(coverage(ind), 's-', 'LineWidth', 1.0 , 'LineStyle', 'none' );
    end
end

ylim([0 1]);
xlim([0 N_examples]);

grid on

ax = gca;
set(ax.XAxis, 'FontSize', 14);
set(ax.YAxis, 'FontSize', 14);


legend('Tube-A', 'Tube-B', 'Tube-C', 'Tube-D', ...
               'SLS-MPC', 'Lumped-Disturbance', 'Offline-Tightening', ...
               'Interpreter', 'latex', 'FontSize', 14, 'Location', 'bestoutside');
xlabel('example number', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('coverage', 'Interpreter', 'Latex', 'FontSize', 20);

mean_coverage = mean(best_coverage_mat)
var_coverage = var(best_coverage_mat)
