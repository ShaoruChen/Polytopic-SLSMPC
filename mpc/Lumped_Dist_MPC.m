classdef Lumped_Dist_MPC < MPC
    %LUMPED_DIST_MPC first lumps both model uncertainty and additive
    %disturbances into a single additive disturbance term. A state
    %feedback/disturbance feedback robust controller is then synthesized
    %based on the lumped disturbance. See 
    %   Monimoy Bujarbaruah, Ugo Rosolia, Yvonne R. Stürz, and Francesco Borrelli. 
    %   "A simple robust MPC for linear systems with parametric and additive uncertainty." 
    %   In 2021 American Control Conference (ACC), pp. 2108-2113. IEEE, 2021.
    %for the details of this methodology. 
    
    properties
        lumped_dist_bound = [];
        % disturbance feedback controller u = M \bar{w} + u_bar
        M = []; u_bar = [];
    end
    
    methods
        function obj = Lumped_Dist_MPC(MPC_args)
            %LUMPED_DIST_MPC Construct an instance of this class
            obj@MPC(MPC_args);
            obj.compute_lumped_dist_bound();
        end
        
        %% compute bounds on the lumped disturbance
        function [bd] = compute_lumped_dist_bound(obj)
            % The lumped disturbance is defined as d_t = Delta_A x_t + Delta_B u_t + w_t
            % This function finds a bound on norm(d_t, inf) given the state
            % and input constraints, the polytopic model uncertainty set and 
            % bounds on w_t. 
            
            % compute inf-norm over-approximation
            norm_type = inf; 
            nx = obj.nx; nu = obj.nu;
            
            % norm bound on x
            X = obj.state_constr;
            if ~isempty(X)
                eye_mat = eye(nx);
                temp = zeros(1, 2*nx);
                for ii = 1:nx
                   e = eye_mat(:, ii);
                   temp(ii) = abs(X.support(e));
                   temp(ii+nx) = abs(-X.support(-e));
                end
                x_bd = max(temp);
            else
               error('Bounded state constraints are required.\n'); 
            end
            
            % norm bound on u
            U = obj.input_constr;
            if ~isempty(U)
                eye_mat = eye(nu);
                temp = zeros(1, 2*nu);
                for ii = 1:nu
                   e = eye_mat(:, ii);
                   temp(ii) = abs(U.support(e));
                   temp(ii+nu) = abs(-U.support(-e));
                end
                u_bd = max(temp);
            else
                error('Bounded input constraints are required.\n'); 
            end
            
            % norm bound on w
            W = obj.system.dist_set;
            if isempty(W)
                w_bd = obj.system.sigma_w;
            else
                eye_mat = eye(nx);
                temp = zeros(1, 2*nx);
                for ii = 1:nx
                   e = eye_mat(:, ii);
                   temp(ii) = abs(W.support(e));
                   temp(ii+nx) = abs(-W.support(-e));
                end
                w_bd = max(temp);
            end
            
            % find uniform bound on the lumped disturbance
            Delta_vertices =  obj.system.Delta_vertices;
            num_vert = length(Delta_vertices);
            
            temp = zeros(1, num_vert);
            for ii = 1:num_vert
                DA = Delta_vertices{ii}.DA; DB = Delta_vertices{ii}.DB;
                temp(ii) = norm(DA, inf)*x_bd + norm(DB, inf)*u_bd + w_bd;
            end
            
            bd = max(temp);
            obj.lumped_dist_bound = bd;
        end
        
        %% solve robust OCP with norm bounded lumped disturbance
        function [sol] = solve(obj, opt, x0)
            %SOLVE solves the robust OCP with norm-bounded lumped
            %disturbances and an LTV state feedback controller.
            
            if nargin < 2
                opt = struct;
                opt.solver = 'mosek';
                opt.verbose = 2;
                opt.problem_type = 'value';
                x0 = obj.x0;
            elseif nargin < 3
                x0 = obj.x0;
            end
            
            if isempty(opt)
                opt = struct;
                opt.solver = 'mosek';
                opt.verbose = 2;
                opt.problem_type = 'value'; 
            end
            
            nx = obj.nx; nu = obj.nu;
            T = obj.horizon;
            
            % get the uniform bound on the lumped disturbance
            if isempty(obj.lumped_dist_bound)
                warning('Bound on the lumped disturbance has to be computed first. \n');
                obj.compute_lumped_dist_bound();
            end            
            sigma_w_unif = obj.lumped_dist_bound;

            uncertain_system = obj.system;
            Q = obj.Q; R = obj.R; Q_T = obj.Q_T;
            
            if isfield(opt, 'problem_type') & strcmp(opt.problem_type, 'optimizer')
                % use a parametric initial condition
                x0 = sdpvar(nx, 1);
            end

            % x = ZA x + ZB u + w           
            Z = kron(diag(ones(1,T),-1), eye(nx)); % block downshift operator
			ZA_block = Z*blkdiag(kron(eye(T), uncertain_system.A), zeros(nx, nx));
            ZB_block = Z*blkdiag(kron(eye(T), uncertain_system.B), zeros(nx, nu));
            
            % parameterize variables
            % an LTV state feedback controller u = M_aug x + ubar_aug is considered
            M = sdpvar( T*nu, T*nx, 'full');
            ubar = sdpvar(T*nu, 1);
            M_aug = [zeros(T*nu, nx) M; zeros(nu, nx) zeros(nu, T*nx)];
            ubar_aug = [ubar; zeros(nu)];
            
            Id = eye(nx*(T+1));
            state_vec_robust = inv(Id- ZA_block)*(Id + ZB_block*M_aug);
            state_vec_robust_x0 = state_vec_robust(:, 1:nx);
            state_vec_robust_w = state_vec_robust(:, nx+1:end);
            state_vec_fixed = inv(Id- ZA_block)*(ZB_block*ubar_aug ) + state_vec_robust_x0*x0;
            
            % state vec = state_vec_fixed + state_vec_robust      
            
            constr = [];
            % structure constraint on M
            for ii = 1:T
               constr = [constr, M((ii-1)*nu+1:ii*nu, (ii-1)*nx+1:end) == zeros(nu, (T-ii+1)*nx)];
            end
            
            % robust state constraints
            Ax = obj.state_constr.A; bx = obj.state_constr.b;
            nFx = size(Ax, 1);
            try
                for ii = 1:T+1
                    for jj = 1:nFx
                    constr = [constr, Ax(jj, :)*state_vec_fixed((ii-1)*nx+1:ii*nx, 1) ...
                                            + norm(Ax(jj, :)*state_vec_robust_w((ii-1)*nx+1:ii*nx,:), 1)*sigma_w_unif <= bx(jj)];
                    end
                end
            catch ME
                sol = struct;
                sol.status = 1; sol.solver_time = nan; sol.solution = [];
                sol.lumped_dist_bound = sigma_w_unif;
                return
            end
            
            % robust terminal constraints
            if ~isempty(obj.terminal_constr)
                Af = obj.terminal_constr.A; bf = obj.terminal_constr.b;
                nFT = size(Af, 1);
                try
                    for jj = 1:nFT
                       constr = [constr, Af(jj, :)*state_vec_fixed(T*nx+1:(T+1)*nx, 1) ...
                                               + norm(Af(jj, :)*state_vec_robust_w(T*nx+1:(T+1)*nx, :), 1)*sigma_w_unif <= bf(jj)];
                    end
                catch
                    sol = struct;
                    sol.status = 1; sol.solver_time = nan; sol.solution = [];
                    sol.lumped_dist_bound = sigma_w_unif;
                    return
                end
            end
            
            % robust input constraint
            Au = obj.input_constr.A; bu = obj.input_constr.b;
            nFu = size(Au, 1);
            
            try
                for ii = 1:T
                    for jj = 1:nFu
                       constr = [constr, Au(jj,:)*ubar((ii-1)*nu+1:ii*nu,1) + norm(Au(jj,:)*M((ii-1)*nu+1:ii*nu, :), 1)*sigma_w_unif <= bu(jj)];
                    end
                end
            catch
                sol = struct;
                sol.status = 1; sol.solver_time = nan; sol.solution = [];
                sol.lumped_dist_bound = sigma_w_unif;
                return
            end
            
            % cost function
            cost = 0;
            Q_mat = kron(eye(T), Q);
            Q_mat = blkdiag(Q_mat, Q_T);
            cost = cost + state_vec_fixed'*Q_mat*state_vec_fixed;
            
            R_mat = kron(eye(T), R);
            cost = cost + ubar'*R_mat*ubar;
            
			ops = sdpsettings('verbose', opt.verbose, 'solver', opt.solver);
      
            if ~isfield(opt, 'problem_type') | strcmp(opt.problem_type, 'value')
                solution = optimize(constr, cost, ops);
            elseif strcmp(opt.problem_type, 'optimizer')
                % the yalmip optimizer is returned
                sol = optimizer(constr, cost, ops, x0, ubar(1:nu));
                return
            end
            
			if solution.problem ~= 0
				warning(['Numerical error detected. Error code ', num2str(solution.problem), '\n']);
                sol = struct;
                sol.solution = solution;
                sol.status = solution.problem;
                sol.solver_time = nan;
                sol.lumped_dist_bound = sigma_w_unif;
                return
			end

			sol = struct;
            sol.solver_time = solution.solvertime;
            sol.solution = solution;
            sol.status = solution.problem;
            
            M_val = value(M);
            M_aug = [zeros(T*nu, nx) M_val; zeros(nu, nx) zeros(nu, T*nx)];
            sol.M = M_aug;
            sol.ubar = value(ubar);
            sol.cost = value(cost);
            sol.x0 = x0;

            sol.lumped_dist_bound = sigma_w_unif;
            
            obj.M = M_aug; obj.u_bar = value(ubar);
            yalmip('clear');
        end
        
        %% simulate open-loop and closed-loop trajectories
        function [traj] = simulate_open_loop_traj(obj, x0, DA_seq, DB_seq, w_seq)
            if isempty(obj.M)
               disp('Solve robust OCP ...');
               obj.solve([], x0);
            end
            
            if nargin < 2
                x0 = obj.x0;
                [DA_seq, DB_seq] = obj.system.sample_model_uncertainty(obj.horizon);
                w_seq = obj.system.sample_disturbance(obj.horizon);
            elseif nargin < 3
                [DA_seq, DB_seq] = obj.system.sample_model_uncertainty(obj.horizon);
                w_seq = obj.system.sample_disturbance(obj.horizon);
            elseif nargin < 5
                w_seq = obj.system.sample_disturbance(obj.horizon);
            end
            
            if isempty(x0)
                x0 = obj.x0;
            end
            
            traj = struct;
            A = obj.system.A; B = obj.system.B;
            
            x_seq = [x0]; u_seq = []; 
            w_bar_seq = []; % w_bar denotes the lumped uncertainty
            
            T = obj.horizon;
            x = x0;
            % feedback controller 
            M = obj.M; u_bar = obj.u_bar;
            nx = obj.nx; nu = obj.nu;
            
            for ii = 1:T
                u = u_bar((ii-1)*nu+1:ii*nu, :);
                if ii > 1
                   u = u + M((ii-1)*nu+1:ii*nu, 1:(ii-1)*nx)*w_bar_seq(:);
                end
                
                x_next = (A + DA_seq{ii})*x + (B + DB_seq{ii})*u + w_seq{ii};
                x_nom = A*x + B*u;
                w_bar = x_next - x_nom;
                w_bar_seq = [w_bar_seq w_bar];
                x_seq = [x_seq x_next];
                u_seq = [u_seq u];
                x = x_next;
            end
            
            % compute trajectory cost
            cost = x_seq(:)'*blkdiag(kron(eye(T), obj.Q), obj.Q_T)*x_seq(:) ... 
                                    + u_seq(:)'*kron(eye(T), obj.R)*u_seq(:);
            
            traj.x = x_seq; traj.u = u_seq; traj.cost = cost;
            traj.w_bar = w_bar_seq;
            traj.DA = DA_seq; traj.DB = DB_seq; traj.w = w_seq;
            traj.method = class(obj);
        end
        
        function [traj_set] = simulate_N_open_loop_traj(obj, N, x0)
            if nargin < 3
                x0 = obj.x0;
            end
            
            traj_set = cell(1, N);
            for ii = 1:N
               traj = obj.simulate_open_loop_traj(x0);
               traj_set{ii} = traj;
            end
        end
        
    end
end

