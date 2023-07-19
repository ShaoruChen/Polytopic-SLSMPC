%% Estimate the feasible domain of the robust MPC methods numerically. 
clear; close all;
data = load('data/MPC_example_data.mat');
MPC_data = data.MPC_data;

%% initial conditions
state_constr = MPC_data.state_constr; 
input_constr = MPC_data.input_constr;
terminal_constr = MPC_data.terminal_constr;
x0_set = terminal_constr.grid(10);

% figure;
% color_1 = [0.8500, 0.3250, 0.0980];
% PlotFcns.plot_polyhedron(state_constr, color_1);
% hold on
% color_2 = [0, 0.4470, 0.7410];
% PlotFcns.plot_polyhedron(terminal_constr, color_2);
% scatter(x0_set(:,1), x0_set(:,2), 'ow', 'filled');

diags_record = struct;

%% SLS MPC
mpc = SLS_MPC(MPC_data);
[diags] = feasibility_evaluation(mpc, x0_set);
diags_record.(diags.method) = diags;

fprintf('%s feasible rate: %.2f, solver time: %.2f', ...
diags.method, diags.feasible_rate, diags.runtime);

figure;
color_1 = [0.8500, 0.3250, 0.0980];
PlotFcns.plot_polyhedron(state_constr, color_1);
hold on
color_2 = [0, 0.4470, 0.7410];
PlotFcns.plot_polyhedron(terminal_constr, color_2);
x0_set = terminal_constr.grid(10);
scatter(x0_set(:,1), x0_set(:,2), 20, 'ow');

feasible_set = diags.feasible_set;
scatter(feasible_set(:,1), feasible_set(:,2), 20, 'oy', 'filled');

xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
title([diags.method, ', T = ', num2str(mpc.horizon)],'Interpreter','none', 'FontSize', 18);

%% tube-MPC-nominal
mpc = Tube_MPC_Nominal(MPC_data);
% mpc.assign_cross_section(nominal_min_dist_inv_set);
mpc.offline_initialization();

[diags] = feasibility_evaluation(mpc, x0_set);
diags_record.(diags.method) = diags;

fprintf('%s feasible rate: %.2f, solver time: %.2f', ...
diags.method, diags.feasible_rate, diags.runtime);

figure;
color_1 = [0.8500, 0.3250, 0.0980];
PlotFcns.plot_polyhedron(state_constr, color_1);
hold on
color_2 = [0, 0.4470, 0.7410];
PlotFcns.plot_polyhedron(terminal_constr, color_2);
x0_set = terminal_constr.grid(10);
scatter(x0_set(:,1), x0_set(:,2), 20, 'ow');

feasible_set = diags.feasible_set;
if ~isempty(feasible_set)
    scatter(feasible_set(:,1), feasible_set(:,2), 20, 'oy', 'filled');
end

xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
title([diags.method, ', T = ', num2str(mpc.horizon)],'Interpreter','none', 'FontSize', 18);


%% tube-MPC-homothetic
mpc = Tube_MPC_Homothetic(MPC_data);
mpc.offline_initialization();

[diags] = feasibility_evaluation(mpc, x0_set);
diags_record.(diags.method) = diags;

fprintf('%s feasible rate: %.2f, solver time: %.2f', ...
diags.method, diags.feasible_rate, diags.runtime);

figure;
color_1 = [0.8500, 0.3250, 0.0980];
PlotFcns.plot_polyhedron(state_constr, color_1);
hold on
color_2 = [0, 0.4470, 0.7410];
PlotFcns.plot_polyhedron(terminal_constr, color_2);
x0_set = terminal_constr.grid(10);
scatter(x0_set(:,1), x0_set(:,2), 20, 'ow');

feasible_set = diags.feasible_set;
if ~isempty(feasible_set)
    scatter(feasible_set(:,1), feasible_set(:,2), 20, 'oy', 'filled');
end
xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
title([diags.method, ', T = ', num2str(mpc.horizon)],'Interpreter','none', 'FontSize', 18);


%% tube-MPC-flexible
mpc = Tube_MPC_Flexible(MPC_data);
mpc.offline_initialization(0.99);

[diags] = feasibility_evaluation(mpc, x0_set);
diags_record.(diags.method) = diags;

fprintf('%s feasible rate: %.2f, solver time: %.2f', ...
diags.method, diags.feasible_rate, diags.runtime);

figure;
color_1 = [0.8500, 0.3250, 0.0980];
PlotFcns.plot_polyhedron(state_constr, color_1);
hold on
color_2 = [0, 0.4470, 0.7410];
PlotFcns.plot_polyhedron(terminal_constr, color_2);
x0_set = terminal_constr.grid(10);
scatter(x0_set(:,1), x0_set(:,2), 20, 'ow');
feasible_set = diags.feasible_set;
if ~isempty(feasible_set)
    scatter(feasible_set(:,1), feasible_set(:,2), 20, 'oy', 'filled');
end

xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
title([diags.method, ', T = ', num2str(mpc.horizon)],'Interpreter','none', 'FontSize', 18);


%% tube MPC
mpc = Tube_MPC(MPC_data);
[nominal_min_dist_inv_set, converge_1] = mpc.initialize_cross_section(20);

mpc.offline_initialization();

[diags] = feasibility_evaluation(mpc, x0_set);
diags_record.(diags.method) = diags;

fprintf('%s feasible rate: %.2f, solver time: %.2f', ...
diags.method, diags.feasible_rate, diags.runtime);

figure;
color_1 = [0.8500, 0.3250, 0.0980];
PlotFcns.plot_polyhedron(state_constr, color_1);
hold on
color_2 = [0, 0.4470, 0.7410];
PlotFcns.plot_polyhedron(terminal_constr, color_2);
x0_set = terminal_constr.grid(10);
scatter(x0_set(:,1), x0_set(:,2), 20, 'ow');

feasible_set = diags.feasible_set;
if ~isempty(feasible_set)
    scatter(feasible_set(:,1), feasible_set(:,2), 20, 'oy', 'filled');
end
xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
title([diags.method, ', T = ', num2str(mpc.horizon)],'Interpreter','none', 'FontSize', 18);

%% lumped disturbance MPC
mpc = Lumped_Dist_MPC(MPC_data);

[diags] = feasibility_evaluation(mpc, x0_set);
diags_record.(diags.method) = diags;

fprintf('%s feasible rate: %.2f, solver time: %.2f', ...
diags.method, diags.feasible_rate, diags.runtime);

figure;
color_1 = [0.8500, 0.3250, 0.0980];
PlotFcns.plot_polyhedron(state_constr, color_1);
hold on
color_2 = [0, 0.4470, 0.7410];
PlotFcns.plot_polyhedron(terminal_constr, color_2);
x0_set = terminal_constr.grid(10);
scatter(x0_set(:,1), x0_set(:,2), 20, 'ow');

feasible_set = diags.feasible_set;
if ~isempty(feasible_set)
    scatter(feasible_set(:,1), feasible_set(:,2), 20, 'oy', 'filled');
end
xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
title([diags.method, ', T = ', num2str(mpc.horizon)],'Interpreter','none', 'FontSize', 18);


%% constr-tightening-MPC
mpc = Constr_Tightening_MPC(MPC_data);
[offline_runtime] = mpc.offline_initialization(3);

[diags] = feasibility_evaluation(mpc, x0_set);
diags_record.(diags.method) = diags;

fprintf('%s feasible rate: %.2f, solver time: %.2f', ...
diags.method, diags.feasible_rate, diags.runtime);

figure;
color_1 = [0.8500, 0.3250, 0.0980];
PlotFcns.plot_polyhedron(state_constr, color_1);
hold on
color_2 = [0, 0.4470, 0.7410];
PlotFcns.plot_polyhedron(terminal_constr, color_2);
x0_set = terminal_constr.grid(10);
scatter(x0_set(:,1), x0_set(:,2), 20, 'ow');

feasible_set = diags.feasible_set;
if ~isempty(feasible_set)
    scatter(feasible_set(:,1), feasible_set(:,2), 20, 'oy', 'filled');
end
xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
title([diags.method, ', T = ', num2str(mpc.horizon)],'Interpreter','none', 'FontSize', 18);

%% summarize success rates
field_names = fieldnames(diags_record);
rate_record = struct;
for ii = 1:length(field_names)
    rate_record.(field_names{ii}) = diags_record.(field_names{ii}).feasible_rate;
end

field_names = fieldnames(diags_record);
runtime_feasible_record = struct;
for ii = 1:length(field_names)
    runtime_feasible_record.(field_names{ii}) = diags_record.(field_names{ii}).avg_runtime_feasible;
end

