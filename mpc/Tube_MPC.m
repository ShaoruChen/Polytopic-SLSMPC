classdef Tube_MPC < MPC
    %TUBE_MPC implemented the RMPC method from 
    %   Wilbur Langson, Ioannis Chryssochoos, S. V. Rakovi?, and David Q. Mayne. 
    %   "Robust model predictive control using tubes." Automatica 40, no. 1 (2004): 125-133.
    % This method can handle both polytopic model uncertainty and additive
    % disturbances. 

    properties
        % cross section: Polyhedron
        cross_section = []; 
        % label denotes the type of the tube cross section
        label = [];
        % controller parameters
        tube = []; u_tube = []; u0 = []; tube_vertices = [];
    end
    
    methods
        function obj = Tube_MPC(MPC_args, tube_args)
            %TUBE_MPC Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                tube_args = [];
            end
            
            obj@MPC(MPC_args);
            
            if ~isempty(tube_args)
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
             
             obj.cross_section = cross_section;
             obj.label = label;
        end
        
        %% initialize cross section
        function [cross_sec, is_converge] = initialize_cross_section(obj, N)
           % Compute the minimal disturbance invariant set for the nominal
           % LTI dynamics as the cross section template.
           
           if nargin < 2
               N = 20;
           end
           
%            obj.system.find_K_LQR(obj.Q, obj.R);
           obj.system.local_stabilizing_controller();
           
%            [cross_sec, is_converge] = obj.system.min_nominal_inv_set([], [], N);
           [cross_sec, is_converge] = obj.system.min_nominal_inv_set_iter(N);

           if ~is_converge
              warning('Minimal disturbance invariant set iterations have not converged yet.');
           else
               disp('Minimal disturbance invariant set found.');
           end
           obj.cross_section = cross_sec;
        end
        
        function [] = offline_initialization(obj)
            if isempty(obj.cross_section)
               obj.initialize_cross_section();
            end
        end
            
        %% solve robust OCP using tubes
        function [sol] = solve(obj, opt, x0, cross_sec)
            if nargin < 2
                opt = struct;
                opt.solver = 'mosek';
                opt.verbose = 2;
                opt.problem_type = 'value';
                x0 = obj.x0;
                cross_sec = obj.cross_section;
            elseif nargin < 3
                x0 = obj.x0;
                cross_sec = obj.cross_section;
            elseif nargin < 4
                cross_sec = obj.cross_section;
            end
               
            if isempty(opt)
                opt = struct;
                opt.solver = 'mosek';
                opt.verbose = 2;
                opt.problem_type = 'value'; 
            end
            
            if isempty(cross_sec)
               cross_sec = obj.cross_section; 
            end
            
            uncertain_system = obj.system;
            Ahat = uncertain_system.A; Bhat = uncertain_system.B;
            nx = obj.nx; nu = obj.nu;
            Q = obj.Q; R = obj.R; Q_T = obj.Q_T;
            T = obj.horizon;
            
            Xc = obj.state_constr; 
            Uc = obj.input_constr;
            Xf = obj.terminal_constr;

            if isfield(opt, 'problem_type') & strcmp(opt.problem_type, 'optimizer')
                % use a parametric initial condition
                x0 = sdpvar(nx, 1);
            end
            
            % process additive disturbances
            if isempty(uncertain_system.dist_set)
                sigma_w = uncertain_system.sigma_w;
                Wc = Polyhedron([eye(nx); -eye(nx)], [sigma_w*ones(nx,1); sigma_w*ones(nx, 1)]);
            else
                Wc = uncertain_system.dist_set;
            end
            
            Aw = Wc.A; bw = Wc.b;

            if isempty(cross_sec)
                if isempty(obj.cross_section)
                    obj.initialize_cross_section();
                    cross_sec = obj.cross_section;
                else
                    cross_sec = obj.cross_section;
                end
            end
            
            cross_sec_vertices = cross_sec.V;
            assert(size(cross_sec_vertices,2) == nx);
            
            nJ = size(cross_sec_vertices, 1);
            % set of vertices of Z
            cross_sec_Vset = mat2cell(cross_sec_vertices', [size(cross_sec_vertices, 2)], ones(1, nJ));
            
            % find the H representation of Z
            Az = cross_sec.A; bz = cross_sec.b;
            assert(isempty(cross_sec.Ae));

            % for computing the Minkowski difference
            nz = size(Az, 1);
            max_W = zeros(nz, 1);
            if sigma_w ~= 0
                opts = optimoptions(@linprog,'Display', 'off');
                for ii = 1:nz
                   [~, fval] = linprog(-Az(ii,:)', Aw, bw,[],[],[],[],opts );
                   max_W(ii) = -fval;
                end
            end
            
            % create optimization variables 
            if T > 1
                alpha = sdpvar(ones(1,T), ones(1,T));
            else
               alpha = {sdpvar(1,1)};
            end
                
            if T > 1
                zc = sdpvar(repmat(nx, 1, T), repmat(1, 1, T));
            else
                zc = {sdpvar(repmat(nx, 1, T), repmat(1, 1, T))};
            end
            
            u0 = sdpvar(nu, 1);
            for ii = 1:T
                % models {u_1, u_2, ..., u_T}
                U{ii} = sdpvar(repmat(nu, 1, nJ), repmat(1, 1, nJ));
                % X{i}{j}: the j-th vertex of the cross section X_i
                % X_i contains state x_i in the prediction
                X{ii} = sdpvar(repmat(nx, 1, nJ), repmat(1, 1, nJ));
            end
            
            for ii = 1:T
                for jj = 1:nJ
                   X{ii}{jj} = zc{ii} + alpha{ii}*cross_sec_Vset{jj};
                end
            end
            
            % construct the optimization problem
            fprintf('Constructing constraints ... \n');
            F = [ismember(u0, Uc)];
            % state and input constraints on tubes
            for ii = 1:T
                F = [F, alpha{ii} >= 0];
                for jj = 1:nJ
                   F = [F, ismember(X{ii}{jj}, Xc), ismember(U{ii}{jj}, Uc)];
                   if ii == T & ~isempty(Xf)
                      F = [F, ismember(X{ii}{jj}, Xf)]; 
                   end
                end
            end
            
            % constraints on states x_2, ..., x_T
            Delta_vertices = obj.system.Delta_vertices;
            num_vert = length(Delta_vertices);
            for ii_delta = 1:num_vert
               DA = Delta_vertices{ii_delta}.DA; 
               DB = Delta_vertices{ii_delta}.DB;
               for ii = 1:T-1
                  for jj = 1:nJ
                    F = [F, Az*(Ahat + DA)*X{ii}{jj} + Az*(Bhat + DB)*U{ii}{jj}...
                         - Az*zc{ii+1} + max_W <= alpha{ii+1}*bz];
                  end          
               end
               % constraint on state x_1
               F = [F, Az*(Ahat + DA)*x0 + Az*(Bhat + DB)*u0 ...
                       - Az*zc{1} + max_W <= alpha{1}*bz];
            end
           
            % construct the cost function
            cost_fcn = @(x,u) x'*Q*x + u'*R*u;
            cost_terminal = @(x) x'*Q_T*x;
                    
            fprintf('Constructing objective function... \n');
            % add stage costs
            cost = cost_fcn(x0, u0);
            for ii = 1:T-1
                for jj = 1:nJ
                    cost = cost + cost_fcn(X{ii}{jj}, U{ii}{jj});
                end
            end
            
            % add terminal cost
            for jj = 1:nJ
               cost = cost + cost_terminal(X{T}{jj}); 
            end
       
            cost = cost/(nJ*T);
            
            % solve the robust OCP
            options = sdpsettings('solver', opt.solver, 'verbose', opt.verbose);
            sol = struct;
            
            fprintf('Tube MPC solver started...\n');
            
            if strcmp(opt.problem_type, 'value')
                solution = optimize(F, cost);
            elseif strcmp(opt.problem_type, 'optimizer')
                % the yalmip optimizer is returned
                sol = optimizer(F, cost, options, x0, U{1});
                return
            end
            
            fprintf('Tube MPC optimization finished. \n');

             % extract the tubes
            tube = cell(1, T);
            for ii = 1:T
                tube{ii} = value(zc{ii}) + value(alpha{ii})*cross_sec;
            end
            
            % extract the vertices
            vertices_set = cell(1, T);
            for ii = 1:T
               vertices = cell(1, nJ);
               for jj = 1:nJ
                  vertices{jj} = value(X{ii}{jj}); 
               end
               vertices_set{ii} = vertices;
            end
            
            % extract the control inputs
            u_tubes = cell(1,T);
            for ii = 1:T
                vertex_input = cell(1, nJ);
                for jj = 1:nJ
                    vertex_input{jj} = value(U{ii}{jj});
                end
               u_tubes{ii} = vertex_input; 
            end
            
            % save the solution
            sol.solution = solution;
            sol.solver_time = solution.solvertime;
            sol.status = solution.problem;
            
            sol.cost = value(cost);
            sol.u0 = value(u0);
            sol.x0 = x0;
            sol.tube = tube;   
            sol.u_tubes = u_tubes;
            sol.vertices_set = vertices_set;
            
            alpha_value = zeros(1, T);
            for ii = 1:T
                alpha_value(ii) = value(alpha{ii});
            end
            sol.alpha = alpha_value;  
            
            obj.tube = tube; obj.u_tube = u_tubes; obj.u0 = value(u0);
            obj.tube_vertices = vertices_set;
            yalmip('clear');            
        end
        
         %% simulate open-loop and closed-loop trajectories
        function [traj] = simulate_open_loop_traj(obj, x0, DA_seq, DB_seq, w_seq)
            if isempty(obj.tube)
               disp('Solve robust OCP ...');
               obj.solve([], x0, []);
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
            
            % feedback controller 
            u_tube = obj.u_tube; u0 = obj.u0;
            tube_vertices = obj.tube_vertices;
            
            x_seq = [x0]; u_seq = [u0]; 
            T = obj.horizon;
            x = x0;
            
            nx = obj.nx; nu = obj.nu;
            
            x = x0; 
            for ii = 1:T
               u = u_seq(:, end);
               x = (A + DA_seq{ii})*x + (B + DB_seq{ii})*u + w_seq{ii};

               x_vert = tube_vertices{ii};
               [coef, status] = find_conv_hull_coeff(x, cell2mat(x_vert)', 1);
               if status ~= 0
                   traj_set = [];
                   error('Tube MPC evaluation is infeasible.');
                   return
               end

               % find tube MPC control inputs
               u_next = cell2mat(u_tube{ii})*coef;
               u_seq = [u_seq u_next];
               x_seq = [x_seq x];
            end
            
            u_seq = u_seq(:, 1:end-1);
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
            for ii = 1:obj.horizon
               poly = tube{ii};
               PlotFcns.plot_polyhedron(poly.projection(proj_dim), color_tube, 'FaceAlpha', 0.3);
            end
            
            xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
            ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
            title([class(obj), ', T = ', num2str(obj.horizon)],'Interpreter','none', 'FontSize', 18)
        end
        
        
    end
end

