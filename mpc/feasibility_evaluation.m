function [diags] = feasibility_evaluation(mpc, x0_set)
%FEASIBILITY_EVALUTION Summary of this function goes here
%   Detailed explanation goes here

num_points = size(x0_set, 1);
feasible_set = []; infeasible_set = [];

opt = struct;
opt.solver = 'mosek'; opt.verbose = 0; opt.problem_type = 'optimizer';
[mpc_optimizer] = mpc.solve(opt);

start_time = tic;
runtime_feasible = 0;
for ii = progress(1:num_points)
    init_x = x0_set(ii, :)';
    inner_start_time = tic;
    [sol_value, errorcode] = mpc_optimizer(init_x);

    if errorcode == 0
        feasible_set = [feasible_set; init_x'];
        inner_runtime = toc(inner_start_time);
        runtime_feasible = runtime_feasible + inner_runtime;
    else
        infeasible_set = [infeasible_set; init_x'];
    end
end
runtime = toc(start_time);

diags = struct;
diags.feasible_set = feasible_set;
diags.infeasible_set = infeasible_set;
diags.feasible_rate = size(feasible_set, 1)/num_points;
diags.runtime = runtime;
if ~isempty(feasible_set)
    diags.avg_runtime_feasible = runtime_feasible/size(feasible_set, 1);
else
    diags.avg_runtime_feasible = nan;
end
diags.x0_set = x0_set;
diags.method = class(mpc);

if isprop(mpc, 'label') 
    diags.label = mpc.label; 
else
    diags.label = [];
end

end

