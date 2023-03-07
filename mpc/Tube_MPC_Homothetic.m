classdef Tube_MPC_Homothetic < MPC
    %TUBE_MPC_HOMOTHETIC implements the robust MPC method proposed in
    %   Matthias Lorenzen, Mark Cannon, and Frank Allgöwer. 
    %   "Robust MPC with recursive model update." 
    %   Automatica 103 (2019): 461-471.
    properties
        K = []; % local stabilizing controller u = Kx
        V = []; nv;
        cross_section = []; % tube cross section given by Vx <= 1
        u_v = []; % controller u(t) = Kx(t) + u_v(t) is synthesized
        
        % label denotes the type of the tube cross section
        label = [];
        
        % off-line parameters for constraint satisfaction
        cx = []; cu = []; cT = [];
        
        % off-line parameters for tube containment 
        w_bar = [];
    end
    
    methods
        function obj = Tube_MPC_Homothetic(MPC_args, tube_args)
            %TUBE_MPC_HOMOTHETIC Construct an instance of this class
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
            obj.compute_constraint_tightening_params();
            obj.compute_tube_containment_params();            
        end
        
        %% initialize local controller 
        function [K] = find_feedback_controller(obj)
           K = obj.system.local_stabilizing_controller();
           obj.K = K;
        end
        
        %% initialize tube cross section
        function [V, converge] = find_cross_section(obj, K, mu, Nstep)
            % 
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
        
        %% compute offline parameters for constraint satisfaction
        function [cx, cu, cT] = compute_constraint_tightening_params(obj)
           if isempty(obj.K)
               warning('Stabilizing feedback controller is not given. \n Finding a feedback controller ...');
               obj.find_feedback_controller();
           end
            
           if isempty(obj.cross_section)
              warning('Tube cross section is not given. \n Initialize cross section ...');
              obj.find_cross_section(); 
           end
           
           P = obj.cross_section; V = obj.V; nv = obj.nv;
           K = obj.K;
           
           % parameters associated with the state constraints
           state_constr = obj.state_constr;
           Ax = state_constr.A; bx = state_constr.b;
           cx = zeros(size(Ax, 1), 1);
           for ii = 1:size(Ax,1)
              cx(ii) = P.support(Ax(ii, :)'); 
           end
           
           % parameters associated with the terminal constraints
           terminal_constr = obj.terminal_constr;
           AT = terminal_constr.A; bT = terminal_constr.b;
           cT = zeros(size(AT, 1), 1);
           for ii = 1:size(AT,1)
              cT(ii) = P.support(AT(ii, :)'); 
           end
           
           % parameters associated with the input constraints
           input_constr = obj.input_constr;
           Au = input_constr.A; bu = input_constr.b;
           Au_bar = Au*K;
           cu = zeros(size(Au_bar, 1), 1);
           for ii = 1:size(Au_bar,1)
              cu(ii) = P.support(Au_bar(ii, :)'); 
           end
           
           obj.cx = cx; obj.cu = cu; obj.cT = cT;
            
        end
        
         %% compute offline parameters for tube containment
        function [w_bar] = compute_tube_containment_params(obj)
            
           if isempty(obj.cross_section)
              warning('Tube cross section is not given. \n Initialize cross section ...');
              obj.find_cross_section(); 
           end
           
           P = obj.cross_section; V = obj.V; nv = obj.nv;
                      
           % compute w_bar
           if isempty(obj.system.dist_set)
               nx = obj.nx;
               sigma_w = obj.system.sigma_w;
               W = Polyhedron([eye(nx); -eye(nx)], [sigma_w*ones(nx,1); sigma_w*ones(nx, 1)]);
           else
               W = obj.system.dist_set; 
           end
           
           w_bar = zeros(nv, 1);
           for ii = 1:nv
              vec = V(ii, :);
              w_bar(ii) = W.support(vec');
           end
           
           obj.w_bar = w_bar;
        end
        
        %% solve the robust OCP
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

            if isempty(obj.K)
                warning('Stabilizing feedback controller is not given. \n Finding a feedback controller ...');
                obj.find_feedback_controller();
            end

            if isempty(obj.V)
                warning('Tube cross section is not given. \n Initialize cross section ...');
                obj.find_cross_section();
            end

            K = obj.K;
            V = obj.cross_section.A; nv = obj.nv;
            assert(norm(V - obj.V) <= 1e-3);
            bv = obj.cross_section.b;

            if isempty(obj.cx) | isempty(obj.cu) | isempty(obj.cT)
                obj.compute_constraint_tightening_params();
            end

            if isempty(obj.w_bar)
                obj.compute_tube_containment_params();
            end

            nx = obj.nx; nu = obj.nu;
            A = obj.system.A; B = obj.system.B;
            T = obj.horizon;

            state_constr = obj.state_constr; input_constr = obj.input_constr;
            terminal_constr = obj.terminal_constr;
            Q = obj.Q; R = obj.R; Q_T = obj.Q_T;

            if isfield(opt, 'problem_type') & strcmp(opt.problem_type, 'optimizer')
                % use a parametric initial condition
                x0 = sdpvar(nx, 1);
            end
            
            cross_sec_vertices = obj.cross_section.V;
            nJ = size(cross_sec_vertices, 1);
            % set of vertices of Z
            cross_sec_Vset = mat2cell(cross_sec_vertices', [size(cross_sec_vertices, 2)], ones(1, nJ));
            
             % extract off-line parameters
            cx = obj.cx; cu = obj.cu; cT = obj.cT; 
            w_bar = obj.w_bar;
            
            % model uncertainty parameters
            Delta_vertices = obj.system.Delta_vertices;
            num_vert = length(Delta_vertices);

            % create optimization variables 
            % cross section is given by zc{t} + alpha{t}*obj.cross_section
            alpha = sdpvar(ones(1,T+1), ones(1,T+1));
            zc = sdpvar(repmat(nx, 1, T+1), repmat(1, 1, T+1));
            
            for ii = 1:T
                % models {u_0, u_1, ..., u_T-1}
                U{ii} = sdpvar(repmat(nu, 1, nJ), repmat(1, 1, nJ));
                % X{i}{j}: the j-th vertex of the cross section X_i
                % X_i contains state x_i in the prediction
                X{ii} = sdpvar(repmat(nx, 1, nJ), repmat(1, 1, nJ));
            end
            
            u_v = sdpvar(repmat(nu, 1, T), repmat(nu, 1, T));
            
            % dual variables
            % uncertainty parameter theta satisfies H_theta*theta <= h_theta.
            H_theta = [-eye(num_vert); eye(num_vert); ones(1, num_vert); -ones(1, num_vert)];
            h_theta = [zeros(num_vert, 1); ones(num_vert, 1); 1; -1];
            
            for ii = 1:T
               Lambda{ii} = sdpvar(repmat(nv, 1, nJ), repmat(size(H_theta, 1), 1, nJ)); 
            end
            
            % use nominal quadratic cost function
            cost = 0;  
            constr = [];
            % nominal state and input variables
            x_bar = sdpvar(repmat(nx, 1, T+1), repmat(1, 1, T+1));
            u_bar = sdpvar(repmat(nu, 1, T), repmat(1, 1, T));
           
            constr = [constr, x_bar{1} == x0];
            for ii = 1:T
               constr = [constr, x_bar{ii+1} == A*x_bar{ii} + B*u_bar{ii}];
               constr = [constr, u_bar{ii} == K*x_bar{ii} + u_v{ii}];
               cost = cost + x_bar{ii}'*Q*x_bar{ii} + u_bar{ii}'*R*u_bar{ii};
            end
            cost = cost + x_bar{T+1}'*Q_T*x_bar{T+1};
            
            % construct constraints
            for ii = 1:T
                for jj = 1:nJ
                   constr = [constr, X{ii}{jj} == zc{ii} + alpha{ii}*cross_sec_Vset{jj}];
                   constr = [constr, U{ii}{jj} == K*X{ii}{jj} + u_v{ii}];
                end
            end
            
            % initial constraints
            constr = [constr, V*(x0 - zc{1}) <= alpha{1}*bv];
            
            % tube containment constraints
            for ii = 1:T
               for jj = 1:nJ
                   constr = [constr, Lambda{ii}{jj} >= 0];
                   
                   D_jj = cell(1, num_vert);
                   for kk = 1:num_vert
                       DA = Delta_vertices{kk}.DA; DB = Delta_vertices{kk}.DB;
                       D_jj{kk} = DA*X{ii}{jj} + DB*U{ii}{jj};
                       constr = [constr, V*D_jj{kk} == Lambda{ii}{jj}*H_theta(:,kk)];
                   end

                   constr = [constr, Lambda{ii}{jj}*h_theta + ... 
                                     V*(A*X{ii}{jj} + B*U{ii}{jj} - zc{ii+1}) + ...
                                                         w_bar <= alpha{ii+1}*bv];
               end 
            end
          
            % constraints tightening
            % state constraint
            Ax = state_constr.A; bx = state_constr.b;
            for ii = 1:T
                constr = [constr, Ax*zc{ii} + alpha{ii}*cx <= bx];
            end
            
            % terminal constraints
            AT = terminal_constr.A; bT = terminal_constr.b;
            constr = [constr, AT*zc{T+1} + alpha{T+1}*cT <= bT];
            
            % input constraints
            Au = input_constr.A; bu = input_constr.b;
            for ii = 1:T
                constr = [constr, Au*K*zc{ii} + Au*u_v{ii} + alpha{ii}*cu <= bu];
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
            
            u_v_value = cell(1, T);
            for ii =1:T
               u_v_value{ii} = value(u_v{ii}); 
            end
            sol.K = K; sol.u_v = cell2mat(u_v_value);
            
            sol.V = V; sol.cross_section = obj.cross_section;
            
            alpha_value = cell(1, T+1); zc_value = cell(1, T+1);
            tube_value = cell(1, T+1);
            for ii = 1:T+1
                alpha_value{ii} = value(alpha{ii});
                zc_value{ii} = value(zc{ii});
                tube_value{ii} = value(zc{ii}) + value(alpha{ii})*obj.cross_section;
            end
            
            sol.alpha = cell2mat(alpha_value); sol.zc = cell2mat(zc_value);
            sol.tube = tube_value;
            
            sol.cost = value(cost);
            
            obj.u_v = cell2mat(u_v_value);
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
            tube = sol.tube;
            for ii = 1:obj.horizon+1
               poly = tube{ii};
               try
                  PlotFcns.plot_polyhedron(poly.projection(proj_dim), color_tube, 'FaceAlpha', 0.3);
               catch 
                  warning('Unable to plot the tube cross section. An empty or degenerate polyhedron may be encountered.');
               end
            end
            
            xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
            ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
            title([class(obj), ', T = ', num2str(obj.horizon)],'Interpreter','none', 'FontSize', 18)
        end
        
    end
end

