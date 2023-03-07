classdef Tube_MPC_Flexible < MPC
    %TUBE_MPC_FLEXIBLE Implements the robust MPC method from
    %   Xiaonan Lu, and Mark Cannon. 
    %   "Robust adaptive tube model predictive control." 
    %   In 2019 American Control Conference (ACC), pp. 3695-3701. IEEE, 2019.
    %A tube MPC approach is applied where the normal vectors of the cross
    %section are fixed. We do not consider uncertainty estimation in this
    %code.
    
    properties
        % hyperparameters
        K = []; % local stabilizing controller u = Kx
        V = []; % cross section parameterized by Vx <= alpha with V fixed
        nv; % V is a matrix of size n_v X n-x.
        u_v = []; % controller u = Kx + u_v
        cross_section = [];
        % label denotes the type of the tube cross section
        label = [];
        
        % off-line parameters for solving OCP
        Hx = []; Hu = []; HT = []; % dual variable matrices for state, input and terminal constraints
        Hv = []; % a 1 X nv cell which contains hat{H}_i for i = 1, ... ,nv. 
        w_bar = []; % deviation caused by the additive disturbance
    end
    
    methods
        function obj = Tube_MPC_Flexible(MPC_args, tube_args)
            %TUBE_MPC_FLEXIBLE Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                tube_args = [];
            end
            
            obj@MPC(MPC_args);
            
            if ~isempty(tube_args)
               if isfield(tube_args, 'K')
                   obj.K = tube_args.K;
               end
               
               if isfield(tube_args, 'V')
                   obj.V = tube_args.V;
                   obj.nv = size(V, 1);
               end
               
               if isfield(tube_args, 'cross_section')
                  obj.cross_section = tube_args.cross_section;
                  obj.V = cross_section.A;
                  obj.nv = size(obj.V, 1);
               end
               
            end
            
        end

        function assign_cross_section(obj, cross_section, label)
            if nargin < 3
                label = [];
            end
            
            % cross section has to be normalized as Vx<=1
            Ac = cross_section.A; bc = cross_section.b;
            Ac = Ac./kron(bc, ones(1, size(Ac, 2)));
            cross_section_normalized = Polyhedron(Ac, ones(size(Ac, 1), 1));
            
            obj.cross_section = cross_section_normalized;
            obj.V = cross_section_normalized.A;
            obj.nv = size(obj.V, 1);
            obj.label = label;
        end
        
        %% offline computation
        function [] = offline_initialization(obj, mu, Nstep)
            if nargin < 2
                mu = 0.99;
                Nstep = 20;
            elseif nargin < 3
                Nstep = 20;
            end
                
            K = obj.find_feedback_controller();
            if isempty(obj.cross_section)
                [V, converge] = obj.find_cross_section(K, mu, Nstep);
            end
            obj.compute_constraint_dual_variables();
            obj.compute_tube_dual_variables();
            obj.disturbance_deviation();
        end
        
        %% initialize local controller 
        function [K] = find_feedback_controller(obj)
           K = obj.system.local_stabilizing_controller();
           obj.K = K;
        end
        
        %% initialize tube cross section
        function [cross_section, converge] = find_cross_section(obj, K, mu, Nstep)            
            if nargin < 2
                K = obj.K;
                mu = 0.99;
                Nstep = 20;
            elseif nargin < 3
                mu = 0.99;
                Nstep = 20;
            elseif nargin < 4
                Nstep = 20;
            end
            
            if isempty(K)
               K = obj.find_feedback_controller(); 
            end
 
            Xinit = obj.state_constr;
            Xinit = Xinit;
            [RIS, diagnostic] = obj.system.contractive_robust_forward_invariant_set(Xinit, K, mu, Nstep);
            
            converge = diagnostic.converge;
            if ~diagnostic.converge
                warning('Contractive robust forward invariant set iterations have not converged.');
            end
            
            % normalize the contractive robust forward invariant set
            Av = RIS.A; bv = RIS.b;
            nv = size(Av, 1);
            assert(all(bv >= 1e-4));
            
            for ii = 1:nv
               Av(ii, :) = Av(ii, :)/bv(ii);
            end
            cross_section = Polyhedron(Av, ones(nv, 1));
            obj.cross_section = cross_section;
            
            V = Av; 
            obj.V = V;  obj.nv = nv;
        end
        
        %% compute offline parameters Hx, Hu, HT
        function [Hx, Hu, HT] = compute_constraint_dual_variables(obj)
            if isempty(obj.V)
                warning('Tube cross section is not given. \n Initialize cross section ...');
                obj.find_cross_section();
            end
            
            if isempty(obj.K)
               warning('Stabilizing feedback controller is not given. \n Finding a feedback controller ...');
               obj.find_feedback_controller();
            end

            V = obj.V; nv = obj.nv;
            K = obj.K;
            
            % compute Hx
            state_constr = obj.state_constr;
            Ax = state_constr.A; bx = state_constr.b;
            Hx = zeros(size(Ax, 1), nv);

            for ii = 1:size(Ax, 1)
               vec = Ax(ii, :); 
               dual_vec = obj.constraint_dual_var_subroutine(vec, V);
               Hx(ii, :) = dual_vec;
            end

            % compute HT
            terminal_constr = obj.terminal_constr;
            AT = terminal_constr.A; bT = terminal_constr.b;
            HT = zeros(size(AT, 1), nv);

            for ii = 1:size(AT, 1)
               vec = AT(ii, :); 
               dual_vec = obj.constraint_dual_var_subroutine(vec, V);
               HT(ii, :) = dual_vec;
            end

            % compute Hu
            input_constr = obj.input_constr;
            Au = input_constr.A; bu = input_constr.b;
            Au_bar = Au*K;
            Hu = zeros(size(Au_bar, 1), nv);

            for ii = 1:size(Au_bar, 1)
               vec = Au_bar(ii, :); 
               dual_vec = obj.constraint_dual_var_subroutine(vec, V);
               Hu(ii, :) = dual_vec;
            end

            obj.Hx = Hx; obj.Hu = Hu; obj.HT = HT;
        end
        
        function [dual_vec] = constraint_dual_var_subroutine(obj, vec, V)
           % vec: 1 X nx, row vector encoding the objective max vec*x
           % V: nv X nx, matrix denoting the input set Vx <= alpha
           % This function finds a feasible dual variable in solving 
           %    max vec*x s.t. Vx <= alpha
           % The problem solved has dual variable mu: 1 X nv
           %    min mu*ones s.t. vec = mu*V, mu >= 0
           
           if nargin < 3
               V = obj.V;
           end
           
           nv = size(V, 1);
           nx = size(vec, 2);
           
           mu = sdpvar(1, nv);
           cost = mu*ones(nv, 1);
           constr = [vec == mu*V, mu >= 0];
           ops = sdpsettings('verbose', 0, 'solver', 'mosek');
           solution = optimize(constr, cost, ops);
           if solution.problem ~= 0
             error('Offline computation of dual variables is infeasible.');
           end
           dual_vec = value(mu);
           yalmip('clear');
        end
        
        %% compute offline parameters Hv for tube containment constraints
        function [Hv] = compute_tube_dual_variables(obj)
            if isempty(obj.V)
                warning('Tube cross section is not given. \n Initialize cross section ...');
                obj.find_cross_section();
            end
            
            V = obj.V; nv = obj.nv;
            Hv = cell(1, nv);
            for ii = 1:nv
               vec = V(ii, :);
               H_i = obj.tube_dual_var_subroutine(vec, V);
               Hv{ii} = H_i;
            end
            obj.Hv = Hv;
        end
        
        function [H_i] = tube_dual_var_subroutine(obj, vec, V)
            % max vec*x1 s.t. V*x0 <= alpha, x1 = dynamics(x0)
            if nargin < 3
                V = obj.V; 
            end
            nv = size(V, 1); 
            nx = obj.nx;
            
            if isempty(obj.K)
                warning('Stabilizing feedback controller is not given. \n Finding a feedback controller ...');
                obj.find_feedback_controller();
            end
            
            K = obj.K;

            Delta_vertices = obj.system.Delta_vertices;
            num_vert = length(Delta_vertices);

            cl_Delta_vertices = cell(1, num_vert+1);
            % nominal closed-loop dynamics
            cl_Delta_vertices{1} = obj.system.A + obj.system.B*K;
            for ii = 1:num_vert
            % uncertain closed-loop dynamics
                cl_Delta_vertices{ii+1} = Delta_vertices{ii}.DA + Delta_vertices{ii}.DB*K; 
            end

            H = sdpvar(num_vert+1, nv);
            ub = sdpvar(1);

            delta_mat = zeros(num_vert+1, nx);
            for ii = 1:num_vert+1
                delta_mat(ii, :) = vec*cl_Delta_vertices{ii}; 
            end

            constr = [delta_mat == H*V];

            theta_vert = [ones(num_vert,1) eye(num_vert)];          
            for ii = 1:num_vert
                constr = [constr, theta_vert(ii, :)*H >= 0 ]; 
                constr = [constr, ub >= theta_vert(ii, :)*H*ones(nv, 1)];
            end

            cost = ub;

            ops = sdpsettings('verbose', 0, 'solver', 'mosek');
            solution = optimize(constr, cost, ops);
            if solution.problem ~= 0
            	error('Offline computation of dual variables is infeasible.');
            end
            H_i = value(H);
            yalmip('clear');
            
        end
        
        %% compute bounds on the deviation caused by the additive disturbances
        function [w_bar] = disturbance_deviation(obj)
            if isempty(obj.V)
                warning('Tube cross section is not given. \n Initialize cross section ...');
                obj.find_cross_section();
            end
            
            V = obj.V; nv = obj.nv;
            nx = obj.nx;
            
            if isempty(obj.system.dist_set)
                sigma_w = obj.system.sigma_w;
                W = Polyhedron([eye(nx); -eye(nx)], [sigma_w*ones(nx,1); sigma_w*ones(nx, 1)]);
            else
                W = obj.system.dist_set;
            end
            
            w_bar = zeros(1, nv);
            for ii = 1:nv
               w_bar(ii) = W.support(V(ii,:)');
            end
            
            obj.w_bar = w_bar;
        end
        
        %% solve the robust OCP
        function [sol] = solve(obj, opt, x0)
            %METHOD1 Summary of this method goes here
            %   Online OCP 
            
            if nargin < 2
                opt = struct;
                opt.solver = 'mosek';
                opt.verbose = 2;
                opt.problem_type = 'value';
                x0 = obj.x0;
            elseif nargin < 3
                x0 = obj.x0;
            end

            nx = obj.nx; nu = obj.nu; 
            A = obj.system.A; B = obj.system.B;
            T = obj.horizon;
            Q = obj.Q; R = obj.R; Q_T = obj.Q_T;

            state_constr = obj.state_constr; input_constr = obj.input_constr;
            terminal_constr = obj.terminal_constr;
            
            if isfield(opt, 'problem_type') & strcmp(opt.problem_type, 'optimizer')
                % use a parametric initial condition
                x0 = sdpvar(nx, 1);
            end
            
            if isempty(obj.K)
                warning('Stabilizing feedback controller is not given. \n Finding a feedback controller ...');
                obj.find_feedback_controller();
            end
            
            if isempty(obj.V)
                warning('Tube cross section is not given. \n Initialize cross section ...');
                obj.find_cross_section();
            end
            
            K = obj.K; V = obj.V; nv = obj.nv;
            
            if isempty(obj.Hx) | isempty(obj.Hu) | isempty(obj.HT)
                obj.compute_constraint_dual_variables();
            end
            
            if isempty(obj.Hv)
                obj.compute_tube_dual_variables();
            end
            
            % extract off-line parameters
            Hx = obj.Hx; Hu = obj.Hu; HT = obj.HT; 
            Hv = obj.Hv; w_bar = obj.w_bar;
            
            % optimization variables
            % need to specify the 'full' option to prevent possible
            % symmetry constraints on alpha
            alpha = sdpvar(nv, T+1, 'full'); % tube is given by X_t = {x | Vx <= alpha(:, t)}
            u_v = sdpvar(nu, T, 'full'); % controller is u(t) = Kx + u_v(:, t)
            
            % use nominal quadratic cost function
            cost = 0;     
            constr = [];
            
            % nominal state and input variables
            x_bar = sdpvar(repmat(nx, 1, T+1), repmat(1, 1, T+1));
            u_bar = sdpvar(repmat(nu, 1, T), repmat(1, 1, T));
           
            constr = [constr, x_bar{1} == x0];
            for ii = 1:T
               constr = [constr, x_bar{ii+1} == A*x_bar{ii} + B*u_bar{ii}];
               constr = [constr, u_bar{ii} == K*x_bar{ii} + u_v(:,ii)];
               cost = cost + x_bar{ii}'*Q*x_bar{ii} + u_bar{ii}'*R*u_bar{ii};
            end
            cost = cost + x_bar{T+1}'*Q_T*x_bar{T+1};

            % initial set constraint
            constr = [constr, V*x0 <= alpha(:, 1)];
            
            % state constraints 
            Ax = state_constr.A; bx = state_constr.b;
            for ii = 1:T+1
               constr = [constr, Hx*alpha(:,ii) <= bx]; 
            end
            
            % terminal constraints
            AT = terminal_constr.A; bT = terminal_constr.b;
            constr = [constr, HT*alpha(:, T+1) <= bT];
            
            % input constraints
            Au = input_constr.A; bu = input_constr.b;
            for ii = 1:T
               constr = [constr, Hu*alpha(:,ii) + Au*u_v(:, ii) <= bu]; 
            end
            
            Delta_vertices = obj.system.Delta_vertices;
            num_vert = length(Delta_vertices);
            aug_DB_vertices_tp = cell(1, num_vert+1);
            aug_DB_vertices_tp{1} = obj.system.B'; 
            for ii = 1:num_vert
                aug_DB_vertices_tp{ii+1} = Delta_vertices{ii}.DB';
            end
            aug_DB_mat = cell2mat(aug_DB_vertices_tp);
            
            % tube containment constraints
            theta_vert = [ones(num_vert, 1) eye(num_vert)];          
            for ii = 1:T
                for jj = 1:nv
                   vec = V(jj, :);
                   vec_cell = mat2cell(kron(ones(1, num_vert+1), vec'), [nx], ones(1, num_vert+1));
                   temp = aug_DB_mat*blkdiag(vec_cell{:});
                   temp = temp';
                   for kk = 1:num_vert
                       % list all vertices of uncertainty parameters
                       theta = theta_vert(kk, :);
                       constr = [constr, theta*temp*u_v(:,ii) + w_bar(jj) + theta*Hv{jj}*alpha(:,ii) <= alpha(jj, ii+1)];
                   end
                   
                end
            end
            
            % set invariance constraint on the terminal state
            % Not necessary for the current experiment
%             for jj = 1:nv
%                 vec = V(jj, :);
%                 for kk = 1:num_vert
%                    theta = theta_vert(kk, :);
%                    constr = [constr, alpha(jj, T+1) >= theta*Hv{jj}*alpha(:, T+1) + w_bar(jj)];
%                 end
%             end
            
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
                sol = optimizer(constr, cost, ops, x0, u_v(:,1));
                return
            end
                  
			if solution.problem ~= 0
				warning(['Numerical error detected. Error code ', num2str(solution.problem), '\n']);
            end
            
            sol = struct;
            sol.status = solution.problem;
            sol.solver_time = solution.solvertime;
            sol.solution = solution;
            
            sol.K = K; sol.u_v = value(u_v);
            sol.alpha = value(alpha);
            sol.V = V;
            sol.x0 = x0;
            
            sol.cost = value(cost);
            
            obj.u_v = value(u_v);
            yalmip('clear');
            
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
            u_v = obj.u_v;
            nx = obj.nx; nu = obj.nu;
            
            for ii = 1:T
                u = K*x_seq(:, end) + u_v(:, ii);
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
        
        function [] = plot_simulated_traj(obj, traj_set, sol, proj_dim)
            if nargin < 4
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
            
            % plot the tube
            color_tube = [0.4660 0.6740 0.1880];
            V = sol.V; alpha = sol.alpha;
            for ii = 1:obj.horizon+1
               poly = Polyhedron(V, alpha(:, ii));
               if poly.volume() > 0
                   PlotFcns.plot_polyhedron(poly.projection(proj_dim), color_tube, 'FaceAlpha', 0.3);
               end
            end
            
            xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
            ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
            title([class(obj), ', T = ', num2str(obj.horizon)],'Interpreter','none', 'FontSize', 18)
        end
        
    end
end

