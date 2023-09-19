clear;
%% initialization
nx = 2;
nu = 1;

x0 = [-3; 4];

N_trials = 120;

diags_list = cell(1, N_trials);
converge_list = zeros(1, N_trials);

for ii = progress(1:N_trials)
eps_A = 0.2; sigma_w = 0.2;  eps_B = 0.1;

A = (rand(nx, nx) - 0.5)*4;
B = (rand(nx, nu) - 0.5)*2;

% scale A
A = A/max(abs(eig(A)))*(0.5 + 2*rand(1));

system_params = struct;
system_params.A = A;
system_params.B = B;

system_params.sigma_w = sigma_w;

% construct model uncertainty set
Delta_vertices = cell(1,4);
Delta_1 = struct; Delta_1.DA = eps_A*[1 0; 0 0]; Delta_1.DB = eps_B*[0; 1];
Delta_2 = struct; Delta_2.DA = eps_A*[1 0; 0 0]; Delta_2.DB = -eps_B*[0; 1];
Delta_3 = struct; Delta_3.DA = -eps_A*[1 0; 0 0]; Delta_3.DB = eps_B*[0; 1];
Delta_4 = struct; Delta_4.DA = -eps_A*[1 0; 0 0]; Delta_4.DB = -eps_B*[0; 1];

Delta_vertices{1} = Delta_1; Delta_vertices{2} = Delta_2; 
Delta_vertices{3} = Delta_3; Delta_vertices{4} = Delta_4;
system_params.Delta_vertices = Delta_vertices;

system = Uncertain_LTI_System(system_params);

%% state and input constraints
Uc_vertices = [-4; 4];
Uc = Polyhedron(Uc_vertices);

E_bar = [1 0; 0 1; -1 0; 0 -1];
e_bar = [8; 8; 8; 8];
Xc = Polyhedron(E_bar, e_bar);

% stage cost weights
Q = 10*eye(nx); R = eye(nu); Q_T = Q;

% uncomment to generate maximum robust positive invariant set
opts = struct;
opts.robust = 1; opts.minVol = 0.1;
[RIS, diagnostic] = system.robust_control_invariant_set(Xc, Uc, 100, opts);
if ~diagnostic.converge
    warning('The search for maximal control invariant set has not converged.');
    continue;
else
    fprintf('Maximal robust control invariant set found!');
end

%% construct SLS MPC problem
MPC_data = struct;
MPC_data.system = system;

horizon = 10; 
MPC_data.horizon = horizon;

MPC_data.Q = Q; MPC_data.R = R; MPC_data.Q_T = Q_T;
MPC_data.state_constr = Xc; 
MPC_data.input_constr = Uc;
% MPC_data.terminal_constr = Xc;
MPC_data.terminal_constr = RIS;

MPC_data.x0 = [3; 0];

file_name = ['data\', 'MPC_example_w_', num2str(sigma_w), '_epsA_', num2str(eps_A), '_epsB_', num2str(eps_B) ];
save(file_name, 'MPC_data');

N_grid = 10; N_iter = 100; early_return = 1;
[diags_record, is_converge] = feasibility_comparison(MPC_data, N_grid, N_iter, early_return);

diags_list{ii} = diags_record;
converge_list(ii) = is_converge;

save data/MPC_feasibility_random_4
end

