classdef Uncertain_LTI_System < handle
    %UNCERTAINLTISYSTEM class of an LTI system with uncertainty
    %   x(k+1) = Ax(k) + Bu(k) + D_A x(k) + D_B u(k) + w(k)
    
    properties(SetAccess = public)
        A; B; % nominal system dynamics matrices
        nx; nu; % state and input dimension

        K = []; % u = Kx is a local linear feedback controller
        OLRIS = []; CLRIS = []; % open-loop and closed-loop (with local controller u = Kx) robust invariant set: Polyhedron
        
        % Delta_vertices{i} is a struct that contains Delta_A and Delta_B,
        % i.e., the i-th vertex of the polytopic model uncertainty set
        Delta_vertices = []; 
        
        % separate vertices of Delta_A and Delta_B 
        DA_vert = []; DB_vert = [];
        
        sigma_w = 0; % the ell_inf norm of additive disturbances w(k)
        % polytopic disturbance set of w(k): Polyhedron
        % If dist_set is non-empty (w(k) belongs to a polytope), sigma_w becomes vacuous. 
        dist_set = [];
        min_box_dist = []; % minimal bounding box of the polyhedral disturbance set
    end
    
    methods
        function obj = Uncertain_LTI_System(params)
            %UNCERTAINLTISYSTEM Construct an instance of this class
            %   initialize with given parameters
            %   params: struct
 
            obj.A = params.A;
            obj.B = params.B;
            
            nx = size(params.A, 1);
            nu = size(params.B, 2);
            obj.nx = nx;
            obj.nu = nu;
            
            if isfield(params, 'sigma_w')
                obj.sigma_w = params.sigma_w;
            end
            
            if isfield(params, 'Delta_vertices')
                obj.Delta_vertices = params.Delta_vertices;
            end
           
            if isfield(params, 'dist_set')
                obj.dist_set = params.dist_set;
            end
            
            obj.find_vertices_of_Delta_AB();
            obj.find_min_box_dist_set();
        end
        
        function [min_box] = find_min_box_dist_set(obj)
            nx = obj.nx;
            
            if isempty(obj.dist_set)
                sigma_w = obj.sigma_w;
                min_box = [-sigma_w*ones(nx, 1) sigma_w*ones(nx, 1)];
                obj.min_box_dist = min_box;
            else
                W = obj.dist_set; 
                % find the minimal bounding box of W
                ub = zeros(nx, 1); lb = zeros(nx, 1);
                eye_mat = eye(nx);
                for ii = 1:nx
                    ub(ii) = W.support(eye_mat(:,ii));
                    lb(ii) = -W.support(-eye_mat(:,ii));
                end
                min_box = [lb ub];
                obj.min_box_dist = min_box;
            end
        end
        
        %% find local controller u = Kx through LQR
        function [K, P] = find_K_LQR(obj, Q, R)
            [K, P] = dlqr(obj.A, obj.B, Q, R);
            obj.K = -K;
        end
        
        %% find minimum disturbance invariant set Z_inv for the nominal dynamics (A, B)
        function [output, is_converge] = min_nominal_inv_set_iter(obj, N, K)
            % N: maximum number of iterations
            if nargin < 3
                K = obj.K;
                if isempty(K)
                    error('The local controller u = Kx has to be computed first.');
                end
            end
            
            if isempty(obj.dist_set)
                sigma_w = obj.sigma_w;
                if obj.sigma_w == 0
                   warning('Minimum disturbance invariant set is empty since the disturbances have zero norm.');
                end
                W0 = Polyhedron([eye(obj.nx); -eye(obj.nx)], obj.sigma_w*ones(2*obj.nx, 1));
            else
                W0 = obj.dist_set;
            end
            
            A = obj.A;
            B = obj.B;
            n = size(A, 1); m = size(B, 2);

            Acl = A + B*K;
            W = W0;
            is_converge = 0;
            volList = zeros(1, N+1);
            volList(1) = W.volume();
            for ii = 1:N
               % find disturbance invariant set using iterations
               W = W + Acl^ii*W0; 
               % call minHRep to reduce redundant inequalities
               W.minHRep();
               volList(ii+1) = W.volume();
               if abs(volList(ii+1) - volList(ii)) <1e-2
                   is_converge = 1;
                   break;
               end
            end
            output = W;
        end
        
        function [Fs, is_converge] = min_nominal_inv_set(obj, epsilon, K, N)
            % method from
            % Rakovic, Sasa V., Eric C. Kerrigan, Konstantinos I. Kouramas, and David Q. Mayne. 
            % "Invariant approximations of the minimal robust positively invariant set." 
            % IEEE Transactions on automatic control 50, no. 3 (2005): 406-410.
            
            if nargin < 2
                epsilon = 0.05; K = obj.K; N = 100;
            elseif nargin < 3
                K = obj.K; N = 100;
            elseif nargin < 4
                N = 100;
            end
            
            if isempty(epsilon)
                epsilon = 0.05;
            end
            
            if isempty(K)
                if isempty(obj.K)
                    error('The local controller u = Kx has to be computed first.');
                else
                    K = obj.K;
                end
            end

            if isempty(obj.dist_set)
                sigma_w = obj.sigma_w;
                if obj.sigma_w == 0
                   warning('Minimum disturbance invariant set is empty since the disturbances have zero norm.');
                end
                W0 = Polyhedron([eye(obj.nx); -eye(obj.nx)], obj.sigma_w*ones(2*obj.nx, 1));
            else
                W0 = obj.dist_set;
            end
            
            A = obj.A;
            B = obj.B;
            nx = size(A, 1); nu = size(B, 2);

            Acl = A + B*K;
            
            count = 0; 
            alpha = 10;
            M_val = 10;
            E = eye(nx);
            while(alpha > epsilon/(epsilon + M_val))
                count = count+1;
                alpha = max(W0.support((Acl^count)'*(W0.A)')./W0.b);
                Ms = zeros(2*nx,1);
                for ii = 1:count
                    Ms = Ms+W0.support([(Acl^(ii-1))' -(Acl^(ii-1))']);
                end
                M_val = max(Ms);
                if count > N
                    is_converge = 0;
                    Fs = [];
                    return; 
                end
            end

            Fs = W0;
            for ii =1:count-1
                Fs = Fs+Acl^ii*W0;
            end
            Fs = (1/(1-alpha))*Fs;
            is_converge = 1;
        end
        
        %% compute the minimal robust forward invariant set
        function [Fs, is_converge] = min_robust_inv_set(obj, epsilon, K, N)
            % compute the minimal robust invariant set using the method
            % from:
            % Kouramas, Konstantinos I., Sasa V. Rakovic, Eric C. Kerrigan, John C. Allwright, and David Q. Mayne. 
            % "On the minimal robust positively invariant set for linear difference inclusions." 
            % In Proceedings of the 44th IEEE Conference on Decision and Control, pp. 2296-2301. IEEE, 2005.
            % Due to the combinatorial nature of this method, the iteration
            % upper bound N cannot be chosen too large. 
            
            if nargin < 2
                epsilon = 0.1; K = obj.K; N = 4;
            elseif nargin < 3
                K = obj.K; N = 4;
            elseif nargin < 4
                N = 4;
            end
               
            if isempty(obj.K)
                error('The local controller u = Kx has to be computed first.');
            else
                K = obj.K;
            end

            if isempty(obj.dist_set)
                sigma_w = obj.sigma_w;
                if obj.sigma_w == 0
                   warning('Minimum disturbance invariant set is empty since the disturbances have zero norm.');
                end
                W0 = Polyhedron([eye(obj.nx); -eye(obj.nx)], obj.sigma_w*ones(2*obj.nx, 1));
            else
                W0 = obj.dist_set;
            end
            nx = obj.nx;

            % normalize the RHS of the H-representation of W0
            Aw = W0.A; bw = W0.b;
            Aw = Aw./kron(bw, ones(1, nx));
            W0 = Polyhedron(Aw, ones(size(Aw, 1), 1));
            
            
            Delta_vertices = obj.Delta_vertices;
            num_vert = length(Delta_vertices);
            A = obj.A; B = obj.B;
            cl_vert_cell = cell(1, num_vert);
            for ii = 1:num_vert
                DA = Delta_vertices{ii}.DA; DB = Delta_vertices{ii}.DB;
                % use transpose of the closed-loop dynamics for ease of
                % computation in the main algorithm
                cl_vert_cell{ii} = ((A + DA) + (B + DB)*K)';
            end
            
            alpha = 10;
            M_val = 10;
            count = 0;
            E_mat = [eye(nx) -eye(nx)];
            
            Rs_params = {{eye(nx)}};
            while alpha > epsilon/(epsilon+M_val)
                count = count + 1;
                % compute alpha
                Rs_cell = obj.matrices_multiplication_permutation(cl_vert_cell, count);
                Rs_params{count+1} = Rs_cell;
                
                alpha = max(W0.support(cell2mat(Rs_cell)*kron(eye(size(Rs_cell, 2)),W0.A')));
                
                Ms = zeros(2*nx,1);
                for jj = 1:2*nx
                    Ms_val = 0;
                    for ii = 0:count-1
%                         Ri_cell = obj.matrices_multiplication_permutation(cl_vert_cell, ii);
                        Ri_cell = Rs_params{ii+1};
                        Ri_mat = cell2mat(Ri_cell);
                        Ms_val = Ms_val + max(W0.support(Ri_mat*kron(eye(size(Ri_cell, 2)), E_mat(:,jj))));
                    end
                    Ms(jj, 1) = Ms_val;
                end
                
                M_val = max(Ms);
                if count > N
                    is_converge = 0;
                    Fs = [];
                    return; 
                end
            end
            
            Fs = W0;
            for ii =1:count-1
                Ri_cell = Rs_params{ii+1};
                Rs = obj.compute_Rs(Ri_cell, W0);
                Fs = Fs + Rs;
            end
            Fs = (1/(1-alpha))*Fs;
            is_converge = 1; 
        end
        
        function [output_cell] = matrices_multiplication_permutation(obj, mat_cell, k)
           % mat_cell: 1 x N cell that contains matrices A_1, ..., A_N
           % output_cell: 1 x N^k cell that contains all possible
           % multiplications of k matrices from the mat_cll:
           % A_i_1*A_i_2*...*A_i_k with i_j \in {1, 2, ..., N}.
            
           N = length(mat_cell);
           nx = size(mat_cell{1}, 1);
           if k > 0
              [temp_cell] = obj.matrices_multiplication_permutation(mat_cell, k-1);
              M = length(temp_cell);
              output_cell = cell(1, N*M);
              for ii = 1:N
                  mul_result = mat_cell{ii}*cell2mat(temp_cell);
                  output_cell((ii-1)*M+1:ii*M) = mat2cell(mul_result, [nx], nx*ones(1, M)); 
              end
           else
               output_cell = {eye(nx)};
           end
        end
        
        function [Rs] = compute_Rs(obj, mat_cell, W)
            N = length(mat_cell);
                        
            % Cannot find the convex hull function for the union of
            % Polyhedron. Transform the data into the polytope class for
            % convex hull operation. 
            R = [];
            for ii = 1:N
                R = [R polytope(mat_cell{ii}*W)];
            end
            Rs = hull(R); 
            Rs = Polyhedron(Rs);
        end
        
        %% compute preset for the nominal linear dynamics (A, B)
        function [PreS] = preset(obj, Xc, Uc)
           % compute the preset of Xc for all possible u \in Uc using
           % projection methods. Nominal linear dynamics (A, B) is
           % considered.
          
            F = Xc.A; f = Xc.b;
            G = Uc.A; g = Uc.b;

            A = obj.A; B = obj.B; 
            n = size(A,2); m = size(B, 2);

            nU = size(G, 1); 
            Fbar = [F*A F*B; zeros(nU, n) G]; fbar = [f; g];
            liftedPolyhedron = Polyhedron(Fbar, fbar);

            dims = 1:n; % project onto the space of system states
            PreS = liftedPolyhedron.projection(dims);
        end

        %% compute robust preset
        function [PreS] = robust_preset(obj, Xc, Uc)
            % compute the preset of Xc for all possible u \in Uc
            if isempty(obj.dist_set)
                W = Polyhedron([eye(obj.nx); -eye(obj.nx)], obj.sigma_w*ones(2*obj.nx, 1));
            else
                W = obj.dist_set;
            end
            
            robustXc = Xc - W;
            F = robustXc.A; f = robustXc.b;
            
            G = Uc.A; g = Uc.b;

            A = obj.A; B = obj.B; 

            n = size(A,2); m = size(B, 2);

            nU = size(G, 1); 
               

            Fbar = []; fbar = [];
            
            % in case of no model uncertainty
            constr = [F*A F*B];
            Fbar = [Fbar; constr]; fbar = [fbar; f];
                  
            Delta_vertices = obj.Delta_vertices;
            for ii = 1:length(Delta_vertices)
                DA = Delta_vertices{ii}.DA; DB = Delta_vertices{ii}.DB;
                 constr = [F*(A + DA) F*(B + DB)];
                Fbar = [Fbar; constr]; fbar = [fbar; f];
            end
            
            newRow = [zeros(nU, n) G];
            Fbar = [Fbar; newRow]; fbar = [fbar; g];
            liftedPolyhedron = Polyhedron(Fbar, fbar);
            tic
            fprintf('Remove redundancy of the lifted polyhedron.\n');
            liftedPolyhedron.minHRep();      
            toc
            dims = 1:n; % project onto the space of system states
            PreS = liftedPolyhedron.projection(dims); 
            tic
            fprintf('Remove redundancy of the preset.\n');
            PreS.minHRep();
            toc
        end
        
        %% preset of a linear autonomous autonomous system
        function [PreS] = preset_autonomous(obj, Xc, sysA)
            % find preset of autonomous system x_+  = sysA*x
            F = Xc.A;  f= Xc.b;
            PreS = Polyhedron(F*sysA, f);
        end
        
        %% preset of autonomous system under model uncertainty
        function [PreS] = robust_preset_autonomous(obj, Xc, K, mute_dist, mu)
            % Enumeration of vertices of the uncertainty set is applied.
            if nargin < 3
                K = obj.K;
                mute_dist = 0;
                mu = 1.0;
            elseif nargin < 4
                mute_dist = 0;
                mu = 1.0;
            elseif nargin < 5
                mu = 1.0;
            end
            
            if ~mute_dist
                if isempty(obj.dist_set)
                    W = Polyhedron([eye(obj.nx); -eye(obj.nx)], obj.sigma_w*ones(2*obj.nx, 1));
                else
                    W = obj.dist_set;
                end
                robustXc = Xc - W;
            else
                robustXc = Xc;
            end
            
            F = robustXc.A; f = robustXc.b;

            A = obj.A; B = obj.B;    
            n = size(A,2); m = size(B, 2);
            Fbar = []; fbar = [];

            Delta_vertices = obj.Delta_vertices;
            for ii = 1:length(Delta_vertices)
               DA = Delta_vertices{ii}.DA; DB = Delta_vertices{ii}.DB;
                constr = [F*((A + DA)+(B + DB)*K)/mu];
                    Fbar = [Fbar; constr]; fbar = [fbar; f];
            end
     
            PreS = Polyhedron(Fbar, fbar);
        end

        %% find the robust control invariant set
        function [RIS, diagnostic] = robust_control_invariant_set(obj, Xinit, Uc, Nstep, options)
            % Xc is the given initial set on states; Uc is the constraint
            % on u; Nstep is the maximum step simulated forward.
            % This function computes the robust control invariant set.
            if nargin < 4
                Nstep = 10;
                options = obj.options;
                options.robust = 0;
                options.minVol = 0.5;
            elseif nargin < 5
                options = obj.options;
                options.robust = 0;
                options.minVol = 0.5;
            end
            
            diagnostic.converge = 0;
            diagnostic.samllSetDetected = 0;
            
            volList = zeros(1, Nstep+1);
            volList(1) = Xinit.volume();
            for ii = 1:Nstep
                fprintf('Robust control invariant set iter %d/%d \n', ii, Nstep);
                fprintf('volume: %f \n', volList(ii));
                
                diagnostic.runningStep = ii;
               
%                 figure; Xinit.plot(); 
%                 title(['step = ', num2str(ii), ' total step = ', num2str(Nstep)]);
%                 pause(0.5);

                if options.robust == 0
                    preS = obj.preset(Xinit, Uc);
                else 
                    preS = obj.robust_preset(Xinit, Uc);
                end
                X_new = and(Xinit, preS);
                X_new.minHRep();
                
                save('RIS_temp', 'X_new');
                
                if X_new == Xinit
                    diagnostic.converge = 1;
                    break;
                end
                Xinit = X_new;
                volList(ii+1) = Xinit.volume();
                
                if Xinit.volume() < options.minVol
                    diagnostic.samllSetDetected = 1;
                    disp('RIS iteration: err: small set encountered.');
                    diagnostic.volList = volList;
                    break;
                end
                
                if abs(volList(ii+1) - volList(ii)) < 1e-2
                    diagnostic.converge = 1;
                    disp('RIS iteration: converged.');
                    diagnostic.volList = volList;
                    break;
                end
            end
            RIS = Xinit;
            obj.OLRIS = RIS;
            diagnostic.volList = volList;
        end
        
        %% Robust invariant set for closed-loop systems
        function [RIS, converge] = robust_forward_invariant_set(obj, Xc, Uc, Nstep, options, K)
            % compute the robust invariant set of the closed-loop system
            % A+BK
            if nargin < 6
                K = obj.K;
            end
            
            if nargin < 5
               options = struct;
               options.plot = 0;
               options.robust = 0;
               options.minVol = 0.5;
            end
            
            F = Xc.A; f = Xc.b;
            G = Uc.A; g = Uc.b;
            
            A = obj.A; B = obj.B; 
            
%             [F_unify, G_unify, ~] = convert_Poly2Mat(Xc, Uc);
%             Xc_new = Polyhedron(F_unify + G_unify*K, ones(size(F_unify, 1), 1));

            Fbar = [F; G*K]; fbar = [f; g];
            Xc_new = Polyhedron(Fbar, fbar);
            
            diagnostic.converge = 0;
            diagnostic.samllSetDetected = 0;
            
            Xinit = Xc_new;
            sysA = A + B*K;
            
            volList = zeros(1, Nstep+1);
            volList(1) = Xinit.volume();
            for ii = 1:Nstep
                fprintf('Robust forward invariant set iter %d/%d \n', ii, Nstep);
                diagnostic.runningStep = ii;

                if options.robust ~= 1
                    PreS = obj.preset_autonomous(Xinit, sysA);
                else 
                    PreS = obj.robust_preset_autonomous(Xinit);
                end
                X_new = and(Xinit, PreS);
                volList(ii+1) = X_new.volume();
                if X_new == Xinit
                    diagnostic.converge = 1;
                    diagnostic.volList = volList;
                    break;
                end
                
                Xinit = X_new;
                if Xinit.volume() < options.minVol
                    diagnostic.samllSetDetected = 1;
                    diagnostic.volList = volList;
                    break;
                end
                
                if abs(volList(ii+1) - volList(ii)) < 1e-2
                    diagnostic.converge = 1;
                    diagnostic.volList = volList;
                    break;
                end
                
            end
            RIS = Xinit;
            obj.CLRIS = RIS;
            diagnostic.volList = volList;
            converge = diagnostic.converge;
        end
        
        %% contractive robust forward invariant set
        function [RIS, diagnostic] = contractive_robust_forward_invariant_set(obj, Xinit, K, mu, Nstep)
            % find the maximal mu-contractive robust forward invariant set
            % for x_+ = ((A + Deleta_A) + (B + Delta_B)*K)*x. No state or
            % input constraints are considered. 
            if nargin < 3
                K = obj.K;
                mu = 0.99;
                Nstep = 20;
            elseif nargin < 4
                mu = 0.99;
                Nstep = 20;
            elseif nargin < 5
                Nstep = 20;
            end

            A = obj.A; B = obj.B; 
            
            diagnostic.converge = 0;
            diagnostic.samllSetDetected = 0;
            
            sysA = A + B*K;
            
            volList = zeros(1, Nstep+1);
            X = Xinit;
            volList(1) = X.volume();
            for ii = 1:Nstep
                fprintf('Contractive robust forward invariant set iter %d/%d \n', ii, Nstep);
                diagnostic.runningStep = ii;
                
                % compute robust preset but not considering the
                % additive disturbances w_t
                mute_dist = 1;
                PreS = obj.robust_preset_autonomous(X, K, mute_dist, mu);
                
                X_new = and(X, PreS);
                volList(ii+1) = X_new.volume();
                if X_new == X
                    diagnostic.converge = 1;
                    diagnostic.volList = volList;
                    break;
                end
                
                X = X_new;
                % exit when volume of the set becomes too small
                if X.volume() < 0.01
                    warning('Small set encountered in contractive robust forward invariant set iterations.');
                    diagnostic.samllSetDetected = 1;
                    diagnostic.converge = 0;
                    diagnostic.volList = volList;
                    break;
                end
                
                if abs(volList(ii+1) - volList(ii)) < 1e-2
                    diagnostic.converge = 1;
                    diagnostic.volList = volList;
                    break;
                end
                
            end
            RIS = X;
            obj.CLRIS = RIS;
            diagnostic.volList = volList;
        end
        
        %% compute local stabilizing controller u = K xi
        function [K_value] = local_stabilizing_controller(obj)
           nx = obj.nx; nu = obj.nu;           
           A = obj.A;
           B = obj.B;
           
           G = sdpvar(nx, nx);
           L = sdpvar(nu, nx);
           
           rho = 0.99;
           constr = [];
           
           Delta_vertices = obj.Delta_vertices;
           num_vert = length(Delta_vertices);
           for ii = 1:num_vert
               vertex = Delta_vertices{ii};
               DA = vertex.DA; DB = vertex.DB;
               A_robust = A + DA;
               B_robust = B + DB;
               mat = [rho*G (A_robust*G + B_robust*L)';
                   A_robust*G+B_robust*L  G];
               constr = [constr, mat >= 0];
           end
          
%            cost = -logdet(G);
%            optimize(constr, cost);
           ops = sdpsettings('verbose', 0);
           solution = optimize(constr, [], ops);
           G_value = value(G);
           L_value = value(L);
           K_value = L_value*inv(G_value);
           obj.K = K_value;
        end
        
        %% find all vertices of Delta_A, Delta_B
        function [DA_cell, DB_cell] = find_vertices_of_Delta_AB(obj)
            nx = obj.nx; nu = obj.nu;
            
            Delta_vertices = obj.Delta_vertices;
            num_vert = length(Delta_vertices);
            
            DA_vert = cell(1, num_vert);
            DB_vert = cell(1, num_vert);
            
            for ii = 1:num_vert
              DA_vert{ii} = Delta_vertices{ii}.DA(:);
              DB_vert{ii} = Delta_vertices{ii}.DB(:);
            end
            
            DA_mat = unique(cell2mat(DA_vert)', 'rows')';
            DB_mat = unique(cell2mat(DB_vert)', 'rows')';

            DA_cell = cell(1, size(DA_mat, 2));
            for ii = 1:size(DA_mat, 2)
              DA_cell{ii} = reshape(DA_mat(:, ii), [nx, nx]);
            end
           
            DB_cell = cell(1, size(DB_mat, 2));
            for ii = 1:size(DB_mat, 2)
               DB_cell{ii} = reshape(DB_mat(:, ii), [nx, nu]); 
            end
           
            obj.DA_vert = DA_cell; obj.DB_vert = DB_cell;
        end
        
        %% sample model uncertainty and disturbances
        function [DA_seq, DB_seq] = sample_model_uncertainty(obj, hor, opt)
            % hor: number of the sampled uncertainty matrices
            % opt: 'varying' if time-varying uncertainty matrices are
            % sampled; 'invariant' if time-invariant uncertainty matrices
            % are sampled. 
            
            if nargin < 3
               opt = 'varying';
            end
            
            nx = obj.nx; nu = obj.nu;

            Delta_vertices = obj.Delta_vertices;
            num_vert = length(Delta_vertices);
            
            DA_vert_cell = cell(num_vert, 1); 
            DB_vert_cell = cell(num_vert, 1);
            for ii = 1:num_vert
                DA_vert_cell{ii} = Delta_vertices{ii}.DA;
                DB_vert_cell{ii} = Delta_vertices{ii}.DB;
            end
            DA_vert_mat = cell2mat(DA_vert_cell); DB_vert_mat = cell2mat(DB_vert_cell);
            
            if strcmp(opt, 'varying')
                DA_seq = cell(1, hor); DB_seq = cell(1, hor);
                for ii = 1:hor
                   % uniformly sample from the convex hull of model uncertainty
                   if num_vert == 1
                       DA = Delta_vertices{1}.DA; DB = Delta_vertices{1}.DB;
                   else
                       % get coefficients of the vertices
                       theta = diff(sort([0  rand(1, num_vert - 1) 1]));
                       DA = kron(theta, eye(nx))*DA_vert_mat;
                       DB = kron(theta, eye(nx))*DB_vert_mat;
                   end
                   DA_seq{ii} = DA; DB_seq{ii} = DB;
                end
            elseif strcmp(opt, 'invariant')
                theta = diff(sort([0  rand(1, num_vert - 1) 1]));
                DA = kron(theta, eye(nx))*DA_vert_mat;
                DB = kron(theta, eye(nx))*DB_vert_mat;
                DA_mat = repmat(DA, 1, hor); DB_mat = repmat(DB, 1, hor);
                DA_seq = mat2cell(DA_mat, [nx], nx*ones(1, hor));
                DB_seq = mat2cell(DB_mat, [nx], nu*ones(1, hor));
            else
               error('Unrecognized sampling option.'); 
            end
            
           
        end
        
        function [w_seq] = sample_disturbance(obj, hor)
            % uniformly sample additive disturbances
            nx = obj.nx;
            w_seq = cell(1, hor);
            
            % use rejection sampling
            if isempty(obj.dist_set)
                for ii = 1:hor
                    w_seq{ii} = (2*rand(nx, 1)-1)*obj.sigma_w;
                end
            else
                min_box = obj.min_box_dist;
                lb = min_box(:, 1); ub = min_box(:, 2);
                
                for ii = 1:hor
                    success = 0;
                    for jj = 1:1000
                       w_sample = lb + rand(nx, 1).*(ub - lb);
                       if W.contains(w_sample)
                            w = w_sample;
                            success = 1;
                            break;
                       end
                    end
                    
                    if success
                        w_seq{ii} = w;
                    else
                        w_seq{ii} = zeros(nx, 1); 
                    end
                end
                
            end
            
        end
     
    end
end

