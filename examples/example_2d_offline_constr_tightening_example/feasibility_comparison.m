function [diags_record, is_converge] = feasibility_comparison(MPC_data, N_grid, N_iter, early_return)
%Compare the feasible domain of robust MPC methods
if nargin < 2
    N_grid = 10;
    N_iter = 20; % max number of iterations for finding tube cross sections
    early_return = 1; % exit if tube cross section iterations have not converged
elseif nargin < 3
    N_iter = 20;
    early_return = 0;
elseif nargin < 4
    early_return = 0;
end

%% initial conditions
state_constr = MPC_data.state_constr; 
input_constr = MPC_data.input_constr;
terminal_constr = MPC_data.terminal_constr;
x0_set = terminal_constr.grid(N_grid);

N = 16; % number of methods tried
diags_record = cell(1, N);
count = 1;

%% find tube cross sections
fprintf('Initialize tube cross sections ...\n');
cross_secs = cell(1, 1);
cross_sec_count = 1;

mpc = Tube_MPC(MPC_data);
[nominal_min_dist_inv_set, converge_1] = mpc.initialize_cross_section(N_iter);

if converge_1
    cross_sec = struct; 
    cross_sec.cross_section = nominal_min_dist_inv_set;
    cross_sec.label = 'nom_min_dist_inv';
    cross_secs{cross_sec_count} = cross_sec;
    cross_sec_count = cross_sec_count + 1;
end

% robust forward invariant set
opts = struct;
opts.robust = 1; opts.minVol = 0.1;
Xc = mpc.state_constr; Uc = mpc.input_constr;
[robust_forward_inv_set, converge_2] = mpc.system.robust_forward_invariant_set(Xc, Uc, N_iter, opts);

if converge_2
    cross_sec = struct; 
    cross_sec.cross_section = robust_forward_inv_set;
    cross_sec.label = 'robust_forward_inv';
    cross_secs{cross_sec_count} = cross_sec;
    cross_sec_count = cross_sec_count + 1;
end

mpc = Tube_MPC_Flexible(MPC_data);
[contractive_forward_inv_set, converge_3]= mpc.find_cross_section([], 0.99, N_iter);

if converge_3
    cross_sec = struct; 
    cross_sec.cross_section = contractive_forward_inv_set;
    cross_sec.label = 'robust_contractive';
    cross_secs{cross_sec_count} = cross_sec;
    cross_sec_count = cross_sec_count + 1;
end

if ~(converge_1 | converge_2 | converge_3)
   is_converge = 0;
   if early_return
      diags_record = [];
      is_converge = 0;
      return
   end
else
    is_converge = 1;
end

num_cross_sec = length(cross_secs);

%% tube-MPC
fprintf('Tube-MPC evaluation ...\n');

for jj = 1:num_cross_sec
mpc = Tube_MPC(MPC_data);
cross_section = cross_secs{jj}.cross_section; label = cross_secs{jj}.label;
mpc.assign_cross_section(cross_section, label);
mpc.offline_initialization();

[diags] = feasibility_evaluation(mpc, x0_set);

diags_record{count} = diags;
count = count + 1;

fprintf('%s %s feasible rate: %.2f, solver time: %.2f \n', ...
diags.method, diags.label, diags.feasible_rate, diags.runtime);

yalmip('clear');
end

%% tube-MPC-nominal
fprintf('Tube-MPC-Nominal evaluation ...\n');

for jj = 1:num_cross_sec
mpc = Tube_MPC_Nominal(MPC_data);
cross_section = cross_secs{jj}.cross_section; label = cross_secs{jj}.label;
mpc.assign_cross_section(cross_section, label);
mpc.offline_initialization();

[diags] = feasibility_evaluation(mpc, x0_set);

diags_record{count} = diags;
count = count + 1;

fprintf('%s %s feasible rate: %.2f, solver time: %.2f \n', ...
diags.method, diags.label, diags.feasible_rate, diags.runtime);
yalmip('clear');

end

%% tube-MPC-homothetic
fprintf('Tube-MPC-Homothetic evaluation ...\n');

for jj = 1:num_cross_sec
mpc = Tube_MPC_Homothetic(MPC_data);
cross_section = cross_secs{jj}.cross_section; label = cross_secs{jj}.label;
mpc.assign_cross_section(cross_section, label);
mpc.offline_initialization();

[diags] = feasibility_evaluation(mpc, x0_set);

diags_record{count} = diags;
count = count + 1;

fprintf('%s %s feasible rate: %.2f, solver time: %.2f \n', ...
diags.method, diags.label, diags.feasible_rate, diags.runtime);
yalmip('clear');

end

%% tube-MPC-Flexible
fprintf('Tube-MPC-Flexible evaluation ...\n');

for jj = 1:num_cross_sec
mpc = Tube_MPC_Flexible(MPC_data);
cross_section = cross_secs{jj}.cross_section; label = cross_secs{jj}.label;
mpc.assign_cross_section(cross_section, label);
mpc.offline_initialization();

[diags] = feasibility_evaluation(mpc, x0_set);

diags_record{count} = diags;
count = count + 1;

fprintf('%s %s feasible rate: %.2f, solver time: %.2f \n', ...
diags.method, diags.label, diags.feasible_rate, diags.runtime);
yalmip('clear');

end


%% SLS MPC with full filter parameterization
fprintf('SLS MPC (full) evaluation ...\n');

mpc = SLS_MPC(MPC_data);
[diags] = feasibility_evaluation(mpc, x0_set);

diags_record{count} = diags;
count = count + 1;

fprintf('%s %s feasible rate: %.2f, solver time: %.2f \n', ...
diags.method, diags.label, diags.feasible_rate, diags.runtime);
yalmip('clear');

%% SLS MPC with diagonal filter parameterization
fprintf('SLS MPC (diagonal) evaluation ...\n');

mpc = SLS_MPC(MPC_data, 0);
[diags] = feasibility_evaluation(mpc, x0_set);

diags_record{count} = diags;
count = count + 1;

fprintf('%s %s feasible rate: %.2f, solver time: %.2f \n', ...
diags.method, diags.label, diags.feasible_rate, diags.runtime);

yalmip('clear');

%% lumped disturbance MPC
fprintf('Lumped disturbance MPC evaluation ...\n');

mpc = Lumped_Dist_MPC(MPC_data);
[diags] = feasibility_evaluation(mpc, x0_set);

diags_record{count} = diags;
count = count + 1;

fprintf('%s %s feasible rate: %.2f, solver time: %.2f \n', ...
diags.method, diags.label, diags.feasible_rate, diags.runtime);
yalmip('clear');

%% constr-tightening-MPC
fprintf('Constraint Tightening MPC evaluation ...\n');

mpc = Constr_Tightening_MPC(MPC_data);
[offline_runtime] = mpc.offline_initialization(3);

[diags] = feasibility_evaluation(mpc, x0_set);

diags_record{count} = diags;
count = count + 1;

fprintf('%s %s feasible rate: %.2f, solver time: %.2f \n', ...
diags.method, diags.label, diags.feasible_rate, diags.runtime);
yalmip('clear');

end

