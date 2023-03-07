classdef Constr_Tightening_MPC < MPC
    %CONSTR_TIGHTENING_MPC implements the robust MPC method proposed in
    %   Monimoy Bujarbaruah, Ugo Rosolia, Yvonne R. Stürz, Xiaojing Zhang, and Francesco Borrelli. 
    %   "Robust MPC for LPV systems via a novel optimization-based constraint tightening." 
    %   Automatica 143 (2022): 110459.
    
    properties
        % offline parameters in Eq.(A.5)
        t0;
        % offline parameters in Eq.(13)
        t1; t2; t3; tw;
        % offline parameters in Eq.(20)
        td1; td2; td3; 
        % offline parameters in Eq.(22)
        tdA; tdB;
        % norm used for constraint tightening
        norm_type = 2; 
        dual_norm_type = 2;
        
        max_norm_DA; max_norm_DB;
        max_norm_w;
        
        % disturbance feedback controller
        M = []; u_bar = [];
        
    end
    
    methods
        function obj = Constr_Tightening_MPC(MPC_args, args)
            if nargin < 2
                args = [];
            end
            
            obj@MPC(MPC_args);
            
            if ~isempty(args)
                obj.norm_type = args.norm_type;
                obj.dual_norm_type = 1/(1 - 1/obj.norm_type);
            end
            
            obj.max_norm_DA_DB();
            obj.max_norm_dist();
        end
        
        function [] = assign_norm_type(obj, norm_type)
           obj.norm_type = norm_type;
           obj.dual_norm_type = 1/(1 - 1/obj.norm_type);
        end
        
        %% compute the maximum norm of DA and DB vertices
        function [max_norm_DA, max_norm_DB] = max_norm_DA_DB(obj)
            % find the maximum norm of DA
            DA_vert = obj.system.DA_vert;
            num_DA_vert = length(DA_vert);
            norm_type = obj.norm_type;
            
            max_norm_DA = zeros(1, num_DA_vert);
            for ii = 1:num_DA_vert
               max_norm_DA(ii) = norm(DA_vert{ii}, norm_type); 
            end
            
            max_norm_DA = max(max_norm_DA);
            
            % find the maximum norm of DA
            DB_vert = obj.system.DB_vert;
            num_DB_vert = length(DB_vert);
            
            max_norm_DB = zeros(1, num_DB_vert);
            for ii = 1:num_DB_vert
               max_norm_DB(ii) = norm(DB_vert{ii}, norm_type); 
            end
            
            max_norm_DB = max(max_norm_DB);
           
            obj.max_norm_DA = max_norm_DA;
            obj.max_norm_DB = max_norm_DB;
        end
        
        %% compute the maximum norm of w_t
        function [max_norm_w] = max_norm_dist(obj)
           if isempty(obj.system.dist_set)
               nx = obj.nx;
               sigma_w = obj.system.sigma_w;
               W = Polyhedron([eye(nx); -eye(nx)], [sigma_w*ones(nx,1); sigma_w*ones(nx, 1)]);
           else
               W = obj.system.dist_set; 
           end
           
           vert = W.V;
           num_vert = size(vert, 2);
           norms = zeros(1, num_vert);
           for ii = 1:num_vert
               norms(ii) = norm(vert(ii, :), obj.norm_type);
           end
           
           max_norm_w = max(norms);
           obj.max_norm_w = max_norm_w;
        end
        
        %% offline computation 
        function [t0] = initialize_t0_tw(obj, var, N)
           % Vertex enumeration is required to bound t1 but scales poorly
           % with the horizon. Here we implement the computationally
           % cheaper alternative proposed by Bujarbaruah et al. (Appendix A.4 of
           % https://arxiv.org/abs/2007.00930).
           % N: truncated horizon with which we apply vertex enumeration. 
           if nargin < 3
               N = 3;
           end
           nx = obj.nx; nu = obj.nu;
           
           state_constr = obj.state_constr; input_constr = obj.input_constr;
           terminal_constr = obj.terminal_constr;
           T = obj.horizon;
           
           Ax = state_constr.A; bx = state_constr.b;
           nAx = size(Ax, 1);
           
           AT = terminal_constr.A; bT = terminal_constr.b;
           nAT = size(AT, 1);
           
           % total number of linear constraints on the states
           num_constr = (T-1)*nAx + nAT;
           
           % bounds corresponding to x_1
           t0 = zeros(num_constr, 1);
           
           % bounds corresponding to x_2, ..., x_{T-1}
           t0(1:nAx) = zeros(nAx, 1);
           for ii = 1:T-2
               t_val = zeros(nAx, 1);
               for jj = 1:nAx
                  vec = Ax(jj,:);
                  if ii <= N
                      [ub] = obj.t0_tw_vertex_enumeration(vec, ii, var);
                  else
                      [ub_1] = obj.t0_tw_vertex_enumeration(vec, N, var);
                      [ub_2] = obj.t0_tw_consistency(vec, N+1, ii, var);
                      ub = ub_1 + ub_2;
                  end
                  t_val(jj) = ub;
               end
               t0(ii*nAx+1:(ii+1)*nAx) = t_val;
           end
           
           % bounds corresponding to x_T
           t_val = zeros(nAT, 1);
           for jj = 1:nAT
              vec = AT(jj,:);
              if ii <= N
                  [ub] = obj.t0_tw_vertex_enumeration(vec, ii, var);
              else
                  [ub_1] = obj.t0_tw_vertex_enumeration(vec, N, var);
                  [ub_2] = obj.t0_tw_consistency(vec, N+1, ii, var);
                  ub = ub_1 + ub_2;
              end
              t_val(jj) = ub;
           end
           t0((T-1)*nAx+1:end) = t_val;
        end
        
        function [ub] = t0_tw_vertex_enumeration(obj, vec, N, var)
            % This function finds 
            % ub = max ||vec*[D_N-A^N ... D_1-A]||_dual_norm_type
            % for all D_i \in P_{A_Delta}^i using vertex enumeration. 
            
            dual_norm_type = obj.dual_norm_type;
            nx = obj.nx; nu = obj.nu;
            
            % nominal dynamics
            A = obj.system.A;
            
            % extract vertices on A
            DA_vert = obj.system.DA_vert;
            num_DA_vert = length(DA_vert);
            
            % there are a total of num_DA_vert^(1+2+...+N) vertices 
            M = sum(1:N);
            base = 1:num_DA_vert;
            temp = repmat(base, [M 1]);
            temp = mat2cell(temp, ones(1, M), num_DA_vert);
            vert_set = combvec(temp{:});
            
            total_num_vert = size(vert_set, 2);
            result = zeros(1, total_num_vert);
            for ii = 1:total_num_vert
               vert_ind = vert_set(:, ii);
               temp_mat = zeros(nx, nx*N);
               for jj = 1:N
                   A_mat = eye(nx);
                   DA_ind = vert_ind(sum(1:jj-1)+1:sum(1:jj));
                   m = length(DA_ind);
                   assert(m == jj);
                   
                   for kk = 1:m
                      A_mat = A_mat*(A+DA_vert{DA_ind(kk)}); 
                   end
                   
                   if strcmp(var, 't0')
                       temp_mat(:, sum(N-jj)*nx+1:sum(N-jj)*nx+nx) = A_mat - A^jj;
                   elseif strcmp(var, 'tw')
                       temp_mat(:, sum(N-jj)*nx+1:sum(N-jj)*nx+nx) = A_mat;
                   else
                       error('Variable not recognized.');
                   end
               end
               result(ii) = norm(vec*temp_mat, dual_norm_type);
            end
            ub = max(result); 
        end
        
        function [ub] = t0_tw_consistency(obj, vec, N, M, var)
            % This function finds 
            % ub = max ||vec*[D_M-A^M ... D_{N}-A^{N}]||_dual_norm_type
            % for all D_i \in P_{A_Delta}^i using the
            % consistency property of induced norms. 
            
            norm_type = obj.norm_type;
            dual_norm_type = obj.dual_norm_type;
            
            max_norm_DA = obj.max_norm_DA;
            vec_dual_norm = norm(vec, dual_norm_type);
            
            A = obj.system.A;
            
            ub = 0;
            for ii = N:M
                if strcmp(var, 't0')
                    ub = ub + (norm(A, norm_type) + max_norm_DA)^ii - norm(A, norm_type)^ii;
                elseif strcmp(var, 'tw')
                    ub = ub + (norm(A, norm_type) + max_norm_DA)^ii;
                else
                   error('Variable not recognized.'); 
                end
            end
            
            ub = vec_dual_norm*ub; 
        end
                
        function [tdA, tdB] = initialize_tdA_tdB(obj)
           T = obj.horizon;
           nx = obj.nx; nu = obj.nu;
           
           norm_type = obj.norm_type; dual_norm_type = obj.dual_norm_type;
           
           Ax = obj.state_constr.A; bx = obj.state_constr.b;
           AT = obj.terminal_constr.A; bT = obj.terminal_constr.b;
           
           nAx = size(Ax, 1); nAT = size(AT, 1);
           tdA = zeros((T-1)*nAx + nAT, 1);
           tdB = zeros((T-1)*nAx + nAT, 1);
           
           DA_vert = obj.system.DA_vert;
           DB_vert = obj.system.DB_vert;
           num_DA_vert = length(DA_vert); num_DB_vert = length(DB_vert);
           
           A = obj.system.A;
           A_pow_cell = cell(1, T);
           for ii = 1:T
              A_pow_cell{ii} = A^(T-ii); 
           end
           A_pow_mat = cell2mat(A_pow_cell);
           
           % compute tdA
           for ii = 1:T-1
               for jj = 1:nAx
                  vec = Ax(jj,:);
                  result = zeros(1, num_DA_vert);
                  for kk = 1:num_DA_vert
                      DA  = DA_vert{kk};
                      diag_mat = mat2cell(repmat(DA, 1, ii), [nx], ones(1, ii)*nx);
                      result(kk) = norm(vec*A_pow_mat(:,(T-ii)*nx+1:end)*blkdiag(diag_mat{:}), dual_norm_type);
                  end
                  ub = max(result);
                  tdA((ii-1)*nAx+jj) = ub;
               end
           end
           
           for jj = 1:nAT
              vec = AT(jj,:);
              result = zeros(1, num_DA_vert);
              for kk = 1:num_DA_vert
                  DA  = DA_vert{kk};
                  diag_mat = mat2cell(repmat(DA, 1, T), [nx], ones(1, T)*nx);
                  result(kk) = norm(vec*A_pow_mat*blkdiag(diag_mat{:}), dual_norm_type);
              end
              ub = max(result);
              tdA((T-1)*nAx+jj) = ub;
           end
           
           % compute tdA
           for ii = 1:T-1
               for jj = 1:nAx
                  vec = Ax(jj,:);
                  result = zeros(1, num_DB_vert);
                  for kk = 1:num_DB_vert
                      DB  = DB_vert{kk};
                      diag_mat = mat2cell(repmat(DB, 1, ii), [nx], ones(1, ii)*nu);
                      result(kk) = norm(vec*A_pow_mat(:,(T-ii)*nx+1:end)*blkdiag(diag_mat{:}), dual_norm_type);
                  end
                  ub = max(result);
                  tdB((ii-1)*nAx+jj) = ub;
               end
           end
           
           for jj = 1:nAT
              vec = AT(jj,:);
              result = zeros(1, num_DA_vert);
              for kk = 1:num_DB_vert
                  DB  = DB_vert{kk};
                  diag_mat = mat2cell(repmat(DB, 1, T), [nx], ones(1, T)*nu);
                  result(kk) = norm(vec*A_pow_mat*blkdiag(diag_mat{:}), dual_norm_type);
              end
              ub = max(result);
              tdB((T-1)*nAx+jj) = ub;
           end
           
           obj.tdA = tdA; obj.tdB = tdB;
        end
        
        function [runtime] = offline_initialization(obj, N)
            if nargin < 2
                N = 3;
            end
            
            start_time = tic;
            
            t0 = obj.initialize_t0_tw('t0', N);
            t1 = obj.max_norm_DA*t0;
            t2 = obj.max_norm_DB*t0;
            t3 = norm(obj.system.B, obj.norm_type)*t0;
            tw = obj.initialize_t0_tw('tw', N);
            
            [tdA, tdB] = obj.initialize_tdA_tdB();
            
            % By Eq.(21):
            td1 = tdA + t1; td2 = tdB + t2; td3 = tdB + t2 + t3;
            
            obj.t0 = t0; obj.t1 = t1; obj.t2 = t2; obj.t3 = t3; obj.tw = tw;
            obj.tdA = tdA; obj.tdB = tdB;
            obj.td1 = td1; obj.td2 = td2; obj.td3 = td3;
            runtime = toc(start_time);
        end

        %% solve the robust OCP online
        function [sol] = solve(obj, opt, x0)
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
            
            if obj.horizon == 1
                [sol] = obj.solve_horizon_one(opt, x0);
            else
                [sol] = obj.solve_large_horizon(opt, x0);
            end
        end
        
        function [sol] = solve_horizon_one(obj, opt, x0)
            % solve_horizon_one is applied when horizon = 1.
            if nargin < 2
                opt = struct;
                opt.solver = 'mosek';
                opt.verbose = 2;
                opt.problem_type = 'value';
                x0 = obj.x0;
            elseif nargin < 3
                x0 = obj.x0;
            end
            
            norm_type = inf; constr_norm_type = 1;
            
            % extract OCP parameters
            uncertain_system = obj.system;
            A = uncertain_system.A; B = uncertain_system.B;
            nx = obj.nx; nu = obj.nu;
            
            % extract hyperplane representation of the disturbance set
            if isempty(obj.system.dist_set)
               nx = obj.nx;
               sigma_w = obj.system.sigma_w;
               W = Polyhedron([eye(nx); -eye(nx)], [sigma_w*ones(nx,1); sigma_w*ones(nx, 1)]);
            else
               W = obj.system.dist_set; 
            end
            
            Aw = W.A; bw = W.b;
            
            x_bar = sdpvar(nx, 1);
            x_next = sdpvar(nx, 1);
            u_bar = sdpvar(nu, 1);
            
            AT = obj.terminal_constr.A; bT = obj.terminal_constr.b;
            Au = obj.input_constr.A; bu = obj.input_constr.b;
            
            Lambda_x = sdpvar(size(AT, 1), size(bw, 1));
            
            constr = [];
            constr = [constr, x_bar == x0];
            constr = [constr, x_next == A*x_bar + B*u_bar];
            % state constraints
            constr = [constr, Lambda_x >= 0, AT == Lambda_x*Aw];
            
            DA_vert = obj.system.DA_vert; DB_vert = obj.system.DB_vert;
            num_DA_vert = length(DA_vert); num_DB_vert = length(DB_vert);
            
            Delta_vertices = obj.system.Delta_vertices;
            num_vert = length(Delta_vertices);
            for ii = 1:num_vert
                DA = Delta_vertices{ii}.DA; DB = Delta_vertices{ii}.DB;
                constr = [constr, AT*((A + DA)*x_bar + (B + DB)*u_bar) + Lambda_x*bw <= bT];
            end
            
            % input constraints
            constr = [constr, Au*u_bar <= bu];
            
            % cost function
            Q_T = obj.Q_T; Q = obj.Q; R = obj.R;
            cost = x_bar'*Q*x_bar + x_next'*Q_T*x_next + u_bar'*R*u_bar;
            
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
                sol = optimizer(constr, cost, ops, x0, u_bar(1:nu));
                return
            end
                  
			if solution.problem ~= 0
				warning(['Numerical error detected. Error code ', num2str(solution.problem)]);
            end
            
            sol = struct;
            sol.status = solution.problem;
            sol.solver_time = solution.solvertime;
            sol.solution = solution;
            
            sol.M = [];
            sol.u_bar = value(u_bar);
            sol.cost = value(cost);
            yalmip('clear');
        end
        
        function [sol] = solve_large_horizon(obj, opt, x0)
            %solve_large_horizon is applied when horizon > 1
            
            if nargin < 2
                opt = struct;
                opt.solver = 'mosek';
                opt.verbose = 2;
                opt.problem_type = 'value';
                x0 = obj.x0;
            elseif nargin < 3
                x0 = obj.x0;
            end
            
            norm_type = inf; constr_norm_type = 1;
            
            % extract OCP parameters
            uncertain_system = obj.system;
            A = uncertain_system.A; B = uncertain_system.B;
            nx = obj.nx; nu = obj.nu;
            sigma_w = uncertain_system.sigma_w;
            T = obj.horizon;
            
            Q = obj.Q; Q_T = obj.Q_T; R = obj.R;
            
            state_constr = obj.state_constr; 
            input_constr = obj.input_constr;
            terminal_constr = obj.terminal_constr;
            
            if isfield(opt, 'problem_type') & strcmp(opt.problem_type, 'optimizer')
                % use a parametric initial condition
                x0 = sdpvar(nx, 1);
            end
            
            Delta_vertices = uncertain_system.Delta_vertices;
            num_vert = length(Delta_vertices);
            
            % construct optimization variables
            x_bar = sdpvar(T*nx, 1); % concatenation of \bar{x}_0, ..., \bar{x}_{T-1}
            x_bar_T = sdpvar(nx, 1); % \bar{x}_T
            
            % disturbance feedback controller is given by u = M_aug*w +
            % u_bar
            u_bar = sdpvar(T*nu, 1); % concatenation of \bar{u}_0, ..., \bar{u}_{T-1}
            M = sdpvar((T-1)*nu, (T-1)*nx, 'full');
            M_aug = [zeros(nu, (T-1)*nx); M];
            M_aug = [M_aug zeros(T*nu, nx)];
            
            % construct block matrices
            A_bar = kron(eye(T), A); B_bar = kron(eye(T), B);
            eye_d = kron(eye(T), eye(nx));
            A_1_bar = eye_d;
            for ii = 1:T-1
               A_1_bar = A_1_bar + kron(diag(ones(1, T-ii), -ii), A^ii);
            end
           
            norm_type = obj.norm_type;
            dual_norm_type = obj.dual_norm_type;
            
            max_norm_w = obj.max_norm_w;
            
            % constraints
            constr = [];
            % structure constraint on M
            for ii = 1:T-2
               constr = [constr, M((ii-1)*nu+1:ii*nu, ii*nx+1:end) == zeros(nu, (T-ii-1)*nx)];
            end
            
            td1 = obj.td1; td2 = obj.td2; td3 = obj.td3; tw = obj.tw;
            Ax = state_constr.A; bx = state_constr.b;
            AT = terminal_constr.A; bT = terminal_constr.b;
            
            Fx = blkdiag(kron(eye(T-1), Ax), AT);
            fx = [kron(ones(T-1,1), bx); bT];

            % extract hyperplane representation of the disturbance set
            if isempty(obj.system.dist_set)
               nx = obj.nx;
               sigma_w = obj.system.sigma_w;
               W = Polyhedron([eye(nx); -eye(nx)], [sigma_w*ones(nx,1); sigma_w*ones(nx, 1)]);
            else
               W = obj.system.dist_set; 
            end
            
            Aw = W.A; bw = W.b;
            Hw = kron(eye(T), Aw); hw = kron(ones(T, 1), bw);
            
            % dual variable for state constraints
            Lambda_x = sdpvar(size(Fx, 1), size(hw, 1), 'full');

            % state constraints
            state_constr_RHS = fx - td1*norm(x_bar, norm_type) ...
                                  - td3*norm(M_aug, norm_type)*max_norm_w ...
                                  - td2*norm(u_bar, norm_type) - tw*max_norm_w;
            constr = [constr, Fx*(A_bar*x_bar + B_bar*u_bar) + Lambda_x*hw <= state_constr_RHS];
            constr = [constr, Lambda_x >= 0];
            constr = [constr, Lambda_x*Hw == Fx*(B_bar*M_aug+(A_1_bar - eye_d)*B_bar*M_aug + eye_d)];
            
            % input constraints
            Au = input_constr.A; bu = input_constr.b;
            Fu = kron(eye(T), Au); fu = kron(ones(T, 1), bu);
            
            % dual variables for input constraints
            Lambda_u = sdpvar(size(Fu, 1), size(hw, 1));
            constr = [constr, Lambda_u >= 0];
            constr = [constr, Fu*M_aug == Lambda_u*Hw];
            constr = [constr, Lambda_u*hw + Fu*u_bar <= fu];
            
            % dynamics constraint
            constr = [constr, x_bar(1:nx) == x0];
            for ii = 1:T-1
               constr = [constr, x_bar(ii*nx+1:(ii+1)*nx) == ...
                         A*x_bar((ii-1)*nx+1:ii*nx) + B*u_bar((ii-1)*nu+1:ii*nu)];
            end
            constr = [constr, x_bar_T == A*x_bar((T-1)*nx+1:T*nx) + B*u_bar((T-1)*nu+1:T*nu)];
            
            % cost function
            cost_weight = blkdiag(kron(eye(T), Q), Q_T, kron(eye(T), R));
            basis = [x_bar; x_bar_T; u_bar];
            cost = basis'*cost_weight*basis;
            
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
                sol = optimizer(constr, cost, ops, x0, u_bar(1:nu));
                return
            end
                  
			if solution.problem ~= 0
				warning(['Numerical error detected. Error code ', num2str(solution.problem)]);
            end
            
            sol = struct;
            sol.status = solution.problem;
            sol.solver_time = solution.solvertime;
            sol.solution = solution;
            
            M_val = value(M);
            M_aug = [zeros(nu, (T-1)*nx); M_val];
            M_aug = [M_aug zeros(T*nu, nx)];
            sol.M = M_aug;
            sol.u_bar = value(u_bar);
            sol.cost = value(cost);
            
            obj.M = M_aug; obj.u_bar = value(u_bar);
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
            w_bar_seq = []; % w_bar equals to the actual disturbance signal
            
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
                
                w_bar = w_seq{ii};
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
        
        function [] = plot_simulated_traj(obj, traj_set, sol, proj_dim)
            if nargin < 3
                sol = [];
                proj_dim = [1 2];
            elseif nargin < 4
                proj_dim = [1 2];
            end
            
            figure;
            color_1 = [0.8500, 0.3250, 0.0980];
            PlotFcns.plot_polyhedron(obj.state_constr.projection(proj_dim), color_1);
            hold on
            if ~isempty(obj.terminal_constr)
                color_2 = [0, 0.4470, 0.7410];
                PlotFcns.plot_polyhedron(obj.terminal_constr.projection(proj_dim), color_2);
            end
            
            x0 = obj.x0;
            scatter(x0(proj_dim(1)), x0(proj_dim(2)), 20, 'ow');
            
            PlotFcns.plot_traj_set(traj_set, proj_dim);

            xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
            ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
            title([class(obj), ', T = ', num2str(obj.horizon)],'Interpreter','none', 'FontSize', 18)
        end
        
        
        
    end
end

