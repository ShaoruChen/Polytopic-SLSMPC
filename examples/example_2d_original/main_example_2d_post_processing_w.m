%% eps_A comparison
data = load('data/example_2d_diags_list_w_hor_3.mat');
diags_list_w = data.diags_list;
w_list = [0.1];

%% compare the best coverage
N_eps = length(w_list);
method_cell = {'Tube_MPC', 'Tube_MPC_Homothetic', 'Tube_MPC_Nominal', 'Tube_MPC_Flexible',  ...
               'SLS_MPC', 'Lumped_Dist_MPC', 'Constr_Tightening_MPC'};
num_method = 7;

best_coverage_mat = zeros(N_eps, num_method);
solver_time_mat = zeros(N_eps, num_method);
for ii = 1:N_eps
    diags_record = diags_list_w{ii};
    N = length(diags_record);
    best_coverage = zeros(1, length(method_cell));
    for jj = 1:N
        diags = diags_record{jj};
        if ~isempty(diags)
            method = diags.method;
            ind = find(strcmp(method_cell, method));

            coverage = diags.feasible_rate;
            if coverage >= best_coverage(ind)
                best_coverage(ind) = coverage;
                solver_time_mat(ii, ind) = diags.avg_runtime_feasible;
            end
        end
    end
    best_coverage_mat(ii, :) = best_coverage;
end
     
figure;
hold on 
for ii = 1:4
    plot(w_list, best_coverage_mat(:,ii), 's-.', 'LineWidth', 1.0 );
end

for ii = 5:7
    plot(w_list, best_coverage_mat(:,ii), 's-', 'LineWidth', 1.0 );
end

ylim([0 1]);
xlabel('$\sigma_w$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('coverage', 'Interpreter', 'Latex', 'FontSize', 18);
grid on
legend('Tube-A', 'Tube-B', 'Tube-C', 'Tube-D', ...
               'SLS-MPC', 'Lumped-Disturbance', 'Offline-Tightening', ...
               'Interpreter', 'latex', 'FontSize', 10, 'Location', 'bestoutside');
set(gca, 'XTick',0.1:0.1:0.7);

%% plot solver time
solver_time = solver_time_mat(1:5, :);
solver_time(find(isnan(solver_time))) = 0;
figure;
hbar = bar(solver_time');
set(gca,'YScale','log')
xlabel('method', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('average solver time', 'Interpreter', 'Latex', 'FontSize', 18);
set(hbar, {'DisplayName'}, {'0.1','0.2', '0.3', '0.4', '0.5'}');
legend();

method_names = {'Tube-A', 'Tube-B', 'Tube-C', 'Tube-D', ...
               'SLS-MPC', 'Lumped-Disturbance', 'Offline-Tightening'};
           
set(gca, 'XTickLabel',method_names, 'XTick',1:numel(method_names));


%% evaluate the choice of tube parameters
tube_cell = {'nom_min_dist_inv', 'robust_forward_inv', 'robust_contractive'};
method_cell = {'Tube_MPC', 'Tube_MPC_Homothetic', 'Tube_MPC_Nominal', 'Tube_MPC_Flexible'};

num_method = length(method_cell); 
num_tube = length(tube_cell);

results = cell(1, num_method);

for ii = 1:num_method
   method_name = method_cell{ii};
   coverage_mat = zeros(N_eps, num_tube);
   for jj = 1:N_eps
       coverage = nan(1, num_tube);
       diags_record = diags_list_w{jj};
       N = length(diags_record);
       for kk = 1:N
           diags = diags_record{kk};
           if ~isempty(diags)
               method = diags.method;
               if strcmp(method_name, method)
                   ind = find(strcmp(tube_cell, diags.label));
                   coverage(ind) = diags.feasible_rate;
               end
           end
       end
       coverage_mat(jj,:) = coverage;
   end
   results{ii} = coverage_mat;
end

methods_legend = {'Tube-A', 'Tube-B', 'Tube-C', 'Tube-D'};
for ii = 1:4
    figure; hold on;
    method_name = methods_legend{ii};
    diags = results{ii};
    for jj = 1:3
        plot(w_list, diags(:,jj),  's-', 'LineWidth', 1.0 );
    end
    grid on
    xlabel('$\sigma_w$', 'Interpreter', 'Latex', 'FontSize', 18);
    ylabel('coverage', 'Interpreter', 'Latex', 'FontSize', 18);
    legend('min-inv', 'robust-forward-inv', 'robust-contractive', ...
                   'Interpreter', 'latex', 'FontSize', 14);
    title(['Coverage of ', method_name])
end

%% evaluate the full parameterization of SLS MPC
SLS_MPC_coverage_mat = zeros(N_eps, 2);
SLS_MPC_solver_time = zeros(N_eps, 2);

for ii = 1:N_eps
    diags_record = diags_list_w{ii};
    N = length(diags_record);
    sls_mpc_coverage = zeros(1, length(method_cell));
    for jj = 1:N
        diags = diags_record{jj};
        if ~isempty(diags)
            method = diags.method;
            if strcmp(method, 'SLS_MPC')
                if strcmp(diags.label, 'full')
                    SLS_MPC_coverage_mat(ii, 1) = diags.feasible_rate;
                    SLS_MPC_solver_time(ii, 1) = diags.avg_runtime_feasible;
                elseif strcmp(diags.label, 'diagonal')
                    SLS_MPC_coverage_mat(ii, 2) = diags.feasible_rate;
                    SLS_MPC_solver_time(ii, 2) = diags.avg_runtime_feasible;
                end
            end
        end
    end
end

SLS_MPC_solver_time(find(isnan(SLS_MPC_solver_time))) = 0;
avg_solver_time_full = mean(SLS_MPC_solver_time(:,1));
avg_solver_time_diagonal = mean(SLS_MPC_solver_time(:,2));

figure;
plot(w_list, SLS_MPC_coverage_mat(:,1), 's-', 'LineWidth', 1.5 );
hold on 
plot(w_list, SLS_MPC_coverage_mat(:, 2), 's-', 'LineWidth', 1.5);
ylim([0 1]);
xlabel('$\sigma_w$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('coverage', 'Interpreter', 'Latex', 'FontSize', 18);
grid on
legend('full', 'diagonal', 'Interpreter', 'Latex', 'FontSize', 18);

