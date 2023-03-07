classdef SLS_MPC < MPC
    %SLS_MPC synthesizes a robust LTV state-feedback controller for systems
    %with both polytopic model uncertainty and additive disturbances. 
    
    properties
        % system responses
        Phi_x = []; Phi_u = []; K = [];
        % diagonal variables of the filter
        sigma_seq = [];
        label = 'full';
    end
    
    methods (Access = public)
        function obj = SLS_MPC(params, full)
            %SLS_MPC Construct an instance of this class
            if nargin < 2
                full = 1;
            end
                
            obj@MPC(params);
            
            if ~full
                % use only a diagonal filter Sigma
                obj.label = 'diagonal'; 
            end
        end
        
        %% solve the robust optimal control problem
        function [sol] = solve(obj, opt, x0)
            if nargin < 2
                opt = struct;
                opt.solver = 'mosek';
                opt.verbose = 2;
                opt.problem_type = 'value';
                % whether use a fully parameterized filter
                opt.full_param_filter = 1;
                x0 = obj.x0;
            elseif nargin < 3
                x0 = obj.x0;
            end
            
            if isempty(opt)
                opt = struct;
                opt.solver = 'mosek';
                opt.verbose = 2;
                opt.problem_type = 'value';
                % whether use a fully parameterized filter
                opt.full_param_filter = 1;
            end
            
            if strcmp(obj.label, 'diagonal')
               opt.full_param_filter = 0; 
            end
            
            if isempty(obj.system.dist_set)
                % Solve SLS MPC with polytopic model uncertainty and norm-bounded
                % additive disturbances.
                [sol] = obj.solve_norm_bounded_dist(opt, x0);
            else
                % Solve SLS MPC with polytopic model uncertainty and
                % polyhedral additive disturbances. 
                [sol] = obj.solve_polyhedral_dist(opt, x0);
            end
        end
        
        %% consider norm-bounded additive disturbance
        function [sol] = solve_norm_bounded_dist(obj, opt, x0)
            % Solve SLS MPC with inf-norm-bounded additive disturbances.
            % opt.problem_type = 'value' if we solve for a specific solution.
            % opt.problem_type = 'optimizer' if we want to construct a
            % parametric Yalmip problem.
            
            if nargin < 2
                opt = struct;
                opt.solver = 'mosek';
                opt.verbose = 2;
                opt.problem_type = 'value';
                % whether use a fully parameterized filter
                opt.full_param_filter = 1;
                x0 = obj.x0;
            elseif nargin < 3
                x0 = obj.x0;
            end
            
            norm_type = inf; constr_norm_type = 1;
            
            % extract OCP parameters
            uncertain_system = obj.system;
            Ahat = uncertain_system.A; Bhat = uncertain_system.B;
            nx = obj.nx; nu = obj.nu;
            sigma_w = uncertain_system.sigma_w;
            T = obj.horizon;
            
            Xc = obj.state_constr; 
            Uc = obj.input_constr;
            Xf = obj.terminal_constr;
            
            if isfield(opt, 'problem_type') & strcmp(opt.problem_type, 'optimizer')
                % use a parametric initial condition
                x0 = sdpvar(nx, 1);
            end
            
            Delta_vertices = uncertain_system.Delta_vertices;
            num_vert = length(Delta_vertices);
            
            % construct system response variables 
			Phi_x = sdpvar( (T + 1) * nx, (T + 1) * nx, 'full');
			Phi_u = sdpvar( (T + 1) * nu, (T + 1) * nx, 'full');
            % extract matrix blocks
            Phi = [Phi_x; Phi_u];
  
            % construct the sigma matrix
            sigma_seq = sdpvar(1, T*nx);
            if ~isfield(opt, 'full_param_filter') | opt.full_param_filter
                % use a fully parameterized filter 
                Sigma_mat = sdpvar((T+1)*nx, (T+1)*nx, 'full');
            else
                % TODO: double check if this parameterization works well
                Sigma_mat = diag([ones(1, nx) sigma_seq]);
            end
        
            % Construct the objective function using nominal cost
			sqrtQ = sqrtm(obj.Q); sqrtR = sqrtm(obj.R); sqrtQ_T = sqrtm(obj.Q_T);
            
            Q_c = mat2cell([repmat(sqrtQ, 1, T) sqrtQ_T], [nx], nx*ones(1, T+1));
            Q_mat = blkdiag(Q_c{:});
            
            R_c = mat2cell([repmat(sqrtR, 1, T) zeros(nu, nu)], [nu], nu*ones(1, T+1));
            R_mat = blkdiag(R_c{:});
            
            cost = norm(Q_mat*Phi_x(:, 1:nx)*x0, 2)^2 + norm(R_mat*Phi_u(:, 1:nx)*x0, 2)^2;
            
            % add constraints to the problem
            constr = [];
            constr = [constr, sigma_seq >= 0];
            
            % Structure constraint of Phi_x and Phi_u
			for k = 1 : T
				constr = [constr, Phi_x( (k - 1)*nx + 1: k*nx, k*nx + 1: end) == zeros(nx, (T + 1 - k)*nx)];
            end

			for k = 1 : T
				constr = [constr, Phi_u( (k - 1)*nu + 1: k*nu, k*nx + 1: end) == zeros(nu, (T + 1 - k)*nx)];
            end
            
            if ~isfield(opt, 'full_param_filter') | opt.full_param_filter
                % block lower-triangular constraint on Sigma
                for k = 1 : T
                    constr = [constr, Sigma_mat( (k - 1)*nx + 1: k*nx, k*nx + 1: end) == zeros(nx, (T + 1 - k)*nx)];
                end

                constr = [constr, Sigma_mat(1:nx, 1:nx) == eye(nx)];
                for k = 1:T
                    constr = [constr, Sigma_mat(k*nx+1:(k+1)*nx, k*nx+1:(k+1)*nx) == diag(sigma_seq((k-1)*nx+1:k*nx))];
                end
            end
            
            % block downshift operator
			Z = kron(diag(ones(1,T),-1), eye(nx));
			ZA_block = Z*blkdiag(kron(eye(T), uncertain_system.A), zeros(nx, nx));
            ZB_block = Z*blkdiag(kron(eye(T), uncertain_system.B), zeros(nx, nu));
			Id = eye((T + 1)*nx);
            % add affine constraint
            constr = [constr, (Id - ZA_block)*Phi_x - ZB_block*Phi_u == Sigma_mat];
                
            % state constraints
			if ~isempty(obj.state_constr)
				Fx = obj.state_constr.A; bx = obj.state_constr.b;
				nFx = size(Fx, 1); nbx = length(bx);
                
                for ii = 1:T
                   for jj = 1: nFx
                      f = Fx(jj,:); b = bx(jj);
                      LHS = f*Phi_x((ii-1)*nx+1:ii*nx,1:nx)*x0;
                      for kk = 1:ii-1
                          LHS = LHS + norm(f*Phi_x((ii-1)*nx+1:ii*nx,kk*nx+1:(kk+1)*nx), constr_norm_type);
                      end
                      constr = [constr, LHS <= b];  
                   end
                end
            else
                warning('State constraints are not given. \n');
            end
            
            % terminal constraint   
            
            if ~isempty(obj.terminal_constr) | ~isempty(obj.state_constr)
            
                if ~isempty(obj.terminal_constr)
                    Ft = obj.terminal_constr.A; bt = obj.terminal_constr.b;
                    nFt = size(Ft, 1); nbt = length(bt);
                else
                    Ft = Fx; bt = bx; nFt = nFx; nbt = nbx;
                end
           
                for jj = 1:nFt
                    f = Ft(jj,:); b = bt(jj);
                    LHS = f*Phi_x(T*nx+1:(T+1)*nx,1:nx)*x0;
                    for kk = 1:T
                       LHS = LHS + norm(f*Phi_x(T*nx+1:(T+1)*nx,kk*nx+1:(kk+1)*nx),constr_norm_type);
                    end
                    constr = [constr, LHS <= b];
                end
            else
                warning('Terminal constraints are not given. \n');
            end
           
            % add input constraint
            if ~isempty(obj.input_constr)
				Fu = obj.input_constr.A; bu = obj.input_constr.b;
				nFu = size(Fu, 1); nbu = length(bu);
                
                for ii = 1:T+1
                   for jj = 1: nFu
                      f = Fu(jj,:); b = bu(jj);
                      LHS = f*Phi_u((ii-1)*nu+1:ii*nu,1:nx)*x0;
                      for kk = 1:ii-1
                          LHS = LHS + norm(f*Phi_u((ii-1)*nu+1:ii*nu,kk*nx+1:(kk+1)*nx),constr_norm_type);
                      end
                      constr = [constr, LHS <= b];  
                   end
                end
            else
                warning('Input constraints are not given. \n');
            end
            
            % signal containment constraint
            eye_mat = eye(nx);
            for ii = 1:T
                for ii_delta = 1:num_vert
                    DA = Delta_vertices{ii_delta}.DA;
                    DB = Delta_vertices{ii_delta}.DB;

                    for kk = 1:nx
                        % row standard basis vector
                        e_vec = eye_mat(kk, :);

                        LHS = norm(e_vec*(DA*Phi_x((ii-1)*nx+1:ii*nx, 1:nx) + ...
                                DB*Phi_u((ii-1)*nu+1:ii*nu, 1:nx) - ...
                                Sigma_mat(ii*nx+1:(ii+1)*nx, 1:nx))*x0, norm_type);

                        for jj = 1:ii-1
                            % use constr_norm_type due to Holder's
                            % inequality
                            LHS = LHS + norm(e_vec*(DA*Phi_x((ii-1)*nx+1:ii*nx, jj*nx+1:(jj+1)*nx) +...
                                DB*Phi_u((ii-1)*nu+1:ii*nu, jj*nx+1:(jj+1)*nx) - ...
                                Sigma_mat(ii*nx+1:(ii+1)*nx, jj*nx+1:(jj+1)*nx)), constr_norm_type);
                        end
                        constr = [constr, LHS + sigma_w <= sigma_seq((ii-1)*nx+kk)];
                    end
                end
            end
            
            % solve the problem
            if isfield(opt, 'solver')
                solver = opt.solver;
            else
                solver = 'mosek';
            end
             
            if isfield(opt, 'verbose')
                verbose = opt.verbose;
            else
                verbose = 2;
            end
             
			ops = sdpsettings('verbose', verbose, 'solver', solver);
            
            if strcmp(opt.problem_type, 'value')
                solution = optimize(constr, cost, ops);
            elseif strcmp(opt.problem_type, 'optimizer')
                % the yalmip optimizer is returned
                sol = optimizer(constr, cost, ops, x0, Phi_u(1:nu, 1:nx));
                return
            end
                  
			if solution.problem ~= 0
				warning(['Numerical error detected. Error code ', num2str(solution.problem), '\n']);
			end

			sol = struct;
			Phi_x_val = value(Phi_x); Phi_u_val = value(Phi_u);
            % set values of redundant variables in Phi_u as zero
            Phi_u_val(find(isnan(Phi_u_val)==1)) = 0;
            cost_val = value(cost);
            sigma_seq_val = value(sigma_seq);
            
            sol.Phi_x = Phi_x_val; sol.Phi_u = Phi_u_val;
            K_val = Phi_u_val*inv(Phi_x_val);
            sol.K = K_val;
            
            sol.cost = cost_val; 
            sol.sigma_seq = sigma_seq_val;
            sol.Sigma = value(Sigma_mat);
            
            sol.status = solution.problem;
            sol.solver_time = solution.solvertime;
            sol.solution = solution;
            
            obj.Phi_x = Phi_x_val; obj.Phi_u = Phi_u_val; obj.sigma_seq = sigma_seq_val;
            obj.K = K_val;
            % clear all yalmip variables
            yalmip('clear');
            
        end
        
        %% consider polyhedral additive disturbance
        function [sol] = solve_polyhedral_dist(obj, opt, x0)
            sol = [];
            error('Not implemented yet.');
        end
        
        %% simulate open-loop and closed-loop trajectories
        function [traj] = simulate_open_loop_traj(obj, x0, DA_seq, DB_seq, w_seq)
            if isempty(obj.K)
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
            T = obj.horizon;
            x = x0;
            % feedback controller 
            K = obj.K;
            nx = obj.nx; nu = obj.nu;
            
            for ii = 1:T
                u = K((ii-1)*nu+1:ii*nu, 1:ii*nx)*x_seq(:);
                x_next = (A + DA_seq{ii})*x + (B + DB_seq{ii})*u + w_seq{ii};
                x_seq = [x_seq x_next];
                u_seq = [u_seq u];
                x = x_next;
            end
            
            % compute trajectory cost
            cost = x_seq(:)'*blkdiag(kron(eye(T), obj.Q), obj.Q_T)*x_seq(:) ... 
                                    + u_seq(:)'*kron(eye(T), obj.R)*u_seq(:);
            
            traj.x = x_seq; traj.u = u_seq; traj.cost = cost;
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

