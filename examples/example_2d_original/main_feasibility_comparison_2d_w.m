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
w_list = [0.1];

N_trials = length(w_list);

diags_list = cell(1, N_trials);
converge_list = zeros(1, N_trials);

for ii = progress(1:N_trials)
    sigma_w = w_list(ii);
    eps_A = 0.1;  eps_B = 0.1;
    
    system_params.sigma_w = sigma_w;
    
    % construct model uncertainty set
    DA_cell = {[0 eps_A; eps_A 0], [0 -eps_A; eps_A 0], [0 eps_A; -eps_A 0], [0 -eps_A; -eps_A 0]};
    DB_cell = {[0; eps_B], [0; -eps_B], [eps_B; 0], [-eps_B; 0]};
    num_DA = length(DA_cell);
    num_DB = length(DB_cell);
    
    Delta_vertices = cell(1, num_DA*num_DB);
    for kk = 1:num_DA
        for jj = 1:num_DB
            Delta = struct;
            Delta.DA = DA_cell{kk};
            Delta.DB = DB_cell{jj};
            Delta_vertices{(kk-1)*num_DA+jj} = Delta; 
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
    Q = 10*eye(nx); R = eye(nu); Q_T = Q;
    
    % uncomment to generate maximum robust positive invariant set
    opts = struct;
    opts.robust = 1; opts.minVol = 0.1;
    [RIS, diagnostic] = system.robust_control_invariant_set(Xc, Uc, 200, opts);
    if ~diagnostic.converge
        warning('The search for maximal control invariant set has not converged.');
        keyboard;
    else
        fprintf('Maximal robust control invariant set found!');
    end
    save('data/RIS', 'RIS');
    
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
    
    N_grid = 10; N_iter = 200; early_return = 0;
    [diags_record, is_converge] = feasibility_comparison(MPC_data, N_grid, N_iter, early_return);
    
    diags_list{ii} = diags_record;
    converge_list(ii) = is_converge;
    
    save data/example_2d_diags_list_w_hor_3.mat

end
