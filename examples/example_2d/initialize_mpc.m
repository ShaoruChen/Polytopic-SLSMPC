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
sigma_w = 0.1; eps_A = 0.4; eps_B = 0.1;

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
[RIS, diagnostic] = system.robust_control_invariant_set(Xc, Uc, 50, opts);
if ~diagnostic.converge
    warning('The search for maximal control invariant set has not converged.');
end
save('data/RIS', 'RIS');

% uncomment to generate the robust positive invariant set w.r.t. a local controller u = Kx
% opts = struct;
% opts.robust = 1; opts.minVol = 0.5; opts.plot = 0;

% [K, P] = uncertain_system.find_K_LQR(Q, R);
% [RIS, diagnostic] = uncertain_system.robustInvariantSetClosedLoop(Xc, Uc, 30, opts);

% RIS_data = load('data/RIS.mat');
% terminal_set = RIS_data.RIS;


%% construct SLS MPC problem

MPC_data = struct;
MPC_data.system = system;

horizon = 3; 
MPC_data.horizon = horizon;

MPC_data.Q = Q; MPC_data.R = R; MPC_data.Q_T = Q_T;
MPC_data.state_constr = Xc; 
MPC_data.input_constr = Uc;
% MPC_data.terminal_constr = Xc;
MPC_data.terminal_constr = RIS;

MPC_data.x0 = [3; 0];

% file_name = ['data/MPC_example_w_', num2str(sigma_w), '_epsA_', num2str(eps_A), '_epsB_', num2str(eps_B) ];
file_name = 'data/MPC_example_data';
save(file_name, 'MPC_data');

