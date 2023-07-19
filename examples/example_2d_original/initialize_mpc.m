clear;
%% initialization

nx = 2;
nu = 1;

A = [1 0.15; 0.1 1];
B = [0.1; 1.1];

x0 = [-3; 4];

system_params = struct;
system_params.A = A;
system_params.B = B;

% consider norm bounded disturbances

% parameters for example_1
% sigma_w = 0.1; eps_A = 0.1; eps_B = 0.05;

% parameters for example_2
sigma_w = 0.1; eps_A = 0.2; eps_B = 0.1;

system_params.sigma_w = sigma_w;

DA_cell = {[0 eps_A; eps_A 0], [0 -eps_A; eps_A 0], [0 eps_A; -eps_A 0], [0 -eps_A; -eps_A 0]};
DB_cell = {[0; eps_B], [0; -eps_B], [eps_B; 0], [-eps_B; 0]};
num_DA = length(DA_cell);
num_DB = length(DB_cell);

Delta_vertices = cell(1, num_DA*num_DB);
for ii = 1:num_DA
    for jj = 1:num_DB
        Delta = struct;
        Delta.DA = DA_cell{ii};
        Delta.DB = DB_cell{jj};
        Delta_vertices{(ii-1)*num_DA+jj} = Delta; 
    end
end

system_params.Delta_vertices = Delta_vertices;

system = Uncertain_LTI_System(system_params);

%% state and input constraints
Uc_vertices = [-4; 4];
Uc = Polyhedron(Uc_vertices);

E_bar = [1 0; 0 1; -1 0; 0 -1];
e_bar = [8; 8; 8; 8];
Xc = Polyhedron(E_bar, e_bar);

% stage cost weights
Q = 10*eye(nx); R = 2*eye(nu); Q_T = Q;

% uncomment to generate maximum robust positive invariant set
opts = struct;
opts.robust = 1; opts.minVol = 0.1;
[RIS, diagnostic] = system.robust_control_invariant_set(Xc, Uc, 50, opts);
if ~diagnostic.converge
    warning('The search for maximal control invariant set has not converged.');
end
save('RIS', 'RIS');

% uncomment to generate the robust positive invariant set w.r.t. a local controller u = Kx
% opts = struct;
% opts.robust = 1; opts.minVol = 0.5; opts.plot = 0;

% [K, P] = uncertain_system.find_K_LQR(Q, R);
% [RIS, diagnostic] = uncertain_system.robustInvariantSetClosedLoop(Xc, Uc, 30, opts);

% RIS_data = load('RIS.mat');
% terminal_set = RIS_data.RIS;


%% construct SLS MPC problem

MPC_data = struct;
MPC_data.system = system;

horizon = 10; 
MPC_data.horizon = horizon;

MPC_data.Q = Q; MPC_data.R = R; MPC_data.Q_T = Q_T;
MPC_data.state_constr = Xc; 
MPC_data.input_constr = Uc;
MPC_data.terminal_constr = Xc;
% MPC_data.terminal_constr = RIS;

MPC_data.x0 = [3; 0];

file_name = 'data/MPC_example_data_original';
save(file_name, 'MPC_data');
