classdef PlotFcns 
    % plotting functions
    methods (Static)
       function plot_polyhedron(P, varargin)
            P_reduced = projectPolytope2Plane(P);
            switch numel(varargin)
                case 1
                    fill(P_reduced.V(:, 1), P_reduced.V(:, 2), varargin{1});
                case 3
                    fill(P_reduced.V(:, 1), P_reduced.V(:, 2), varargin{1}, varargin{2}, varargin{3});
                case 5
                    fill(P_reduced.V(:, 1), P_reduced.V(:, 2), varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});

            end
            hold on;
       end
       
       function plot_traj_set(traj_set, proj_dim)
            if nargin < 2
               proj_dim = [1 2];
            end

            N = length(traj_set);
            for ii = 1:N
                plot(traj_set{ii}.x(proj_dim(1),:), traj_set{ii}.x(proj_dim(2),:), 's-', 'LineWidth', 1.5);
                hold on
            end  
       end
       
       function [] = plot_SLS_MPC_traj(sls_mpc, trajSet, title_str, proj_dim)
            %PLOTTRAJ plot the trajectories given in the set trajSet
            if nargin < 3
                title_str = ''
                proj_dim = [1 2];
            elseif nargin < 4
                proj_dim = [1 2];
            end

            T = sls_mpc.horizon;

            figure;
            color_1 = [0.8500, 0.3250, 0.0980];
            % color_1 = [0, 0.4470, 0.7410];
            PlotFcns.plot_polyhedron(sls_mpc.state_constr.projection(proj_dim), color_1, 'LineStyle', 'none');
            hold on
            color_2 = [0, 0.4470, 0.7410];
            % color_2 =  [0.4660 0.6740 0.1880];
            if ~isempty(sls_mpc.terminal_constr)
                PlotFcns.plot_polyhedron(sls_mpc.terminal_constr.projection(proj_dim), color_2, 'LineStyle', 'none');
            end

            N = length(trajSet);
            for ii = 1:N
                plot(trajSet{ii}.x_seq(proj_dim(1),:), trajSet{ii}.x_seq(proj_dim(2),:), 's-', 'LineWidth', 1.5);
                hold on
            end
            x0 = sls_mpc.x0;
            scatter(x0(proj_dim(1)), x0(proj_dim(2)),'oy','filled');

            xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
            ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
            title([title_str, ' trajectories, horizon = ', num2str(T)], 'Interpreter', 'Latex', 'FontSize', 18);

            figure;
            for ii = 1:N
                plot(trajSet{ii}.u_seq, 's-', 'LineWidth', 1.5);
                hold on
            end
            xlabel('$time$', 'Interpreter', 'Latex', 'FontSize', 18);
            ylabel('$u$', 'Interpreter', 'Latex', 'FontSize', 18);
            title([title_str, ' control inputs, horizon = ', num2str(T)], 'Interpreter', 'Latex', 'FontSize', 18);
        end

        function [] = plot_tube_MPC_traj(sls_mpc, tube_sol, traj_set, title_str, proj_dim)
            if nargin < 4
                title_str = ''
                proj_dim = [1 2];
            elseif nargin < 5
                proj_dim = [1 2];
            end
            
            T = sls_mpc.horizon;
            
            figure;
            color_1 = [0.8500, 0.3250, 0.0980];
            PlotFcns.plot_polyhedron(sls_mpc.state_constr.projection(proj_dim), color_1, 'LineStyle', 'none');
            hold on
            color_2 = [0, 0.4470, 0.7410];
            if ~isempty(sls_mpc.terminal_constr)
                PlotFcns.plot_polyhedron(sls_mpc.terminal_constr.projection(proj_dim), color_2, 'LineStyle', 'none');
            end

            for ii = 1:sls_mpc.horizon
                color_1 = [178 102 255]./255;
                PlotFcns.plot_polyhedron(tube_sol.tubes{ii}.projection(proj_dim), color_1, 'FaceAlpha', 0.3);
                hold on
            end

            N = length(traj_set);
            for ii = 1:N
                plot(traj_set{ii}.xi_seq(proj_dim(1),:), traj_set{ii}.xi_seq(proj_dim(2),:), 's-', 'LineWidth', 1.5);
                hold on
            end
            x0 = sls_mpc.x0;
            scatter(x0(proj_dim(1)), x0(proj_dim(2)),'oy','filled');

            xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
            ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
            title([title_str,' trajectories, horizon = ', num2str(T)], 'Interpreter', 'Latex', 'FontSize', 18);

            figure;
            for ii = 1:N
                plot(traj_set{ii}.u_seq, 's-', 'LineWidth', 1.5);
                hold on
            end
            xlabel('$time$', 'Interpreter', 'Latex', 'FontSize', 18);
            ylabel('$u$', 'Interpreter', 'Latex', 'FontSize', 18);
            title([title_str, ' control inputs, horizon = ',num2str(T)], 'Interpreter', 'Latex', 'FontSize', 18);

        end
       
        function [] = plot_uniform_dist_feedback_MPC_traj(sls_mpc, traj_set, title_str, proj_dim)
            if nargin < 3
                title_str = ''
                proj_dim = [1 2];
            elseif nargin < 4
                proj_dim = [1 2];
            end
            
            T = sls_mpc.horizon;
            
            figure;
            color_1 = [0.8500, 0.3250, 0.0980];
            PlotFcns.plot_polyhedron(sls_mpc.state_constr.projection(proj_dim), color_1, 'LineStyle', 'none');
            hold on
            color_2 = [0, 0.4470, 0.7410];
            if ~isempty(sls_mpc.terminal_constr)
                PlotFcns.plot_polyhedron(sls_mpc.terminal_constr.projection(proj_dim), color_2, 'LineStyle', 'none');
            end

            N = length(traj_set);
            for ii = 1:N
                plot(traj_set{ii}.x_seq(proj_dim(1),:), traj_set{ii}.x_seq(proj_dim(2),:), 's-', 'LineWidth', 1.5);
                hold on
            end
            x0 = sls_mpc.x0;
            scatter(x0(proj_dim(1)), x0(proj_dim(2)),'oy','filled');

            xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 18);
            ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 18);
            title([title_str, ' trajectories, horizon = ', num2str(T)], 'Interpreter', 'Latex', 'FontSize', 18);

            figure;
            for ii = 1:N
                plot(traj_set{ii}.u_seq, 's-', 'LineWidth', 1.5);
                hold on
            end
            xlabel('$time$', 'Interpreter', 'Latex', 'FontSize', 18);
            ylabel('$u$', 'Interpreter', 'Latex', 'FontSize', 18);
            title([title_str, ' control inputs, horizon = ', num2str(T)], 'Interpreter', 'Latex', 'FontSize', 18);
        end
    end
end

function P_projected = projectPolytope2Plane(P, dim)
    if nargin < 2
        dim = [1 2];
    end
    
    vert = P.V;
    x_vert = vert(:, dim(1));
    y_vert = vert(:, dim(2));
%     x_vert = round(vert(:, dim(1)), 5);
%     y_vert = round(vert(:, dim(2)), 5);
    idx = convhull(x_vert, y_vert);
    P_projected = Polyhedron([x_vert(idx), y_vert(idx)]);
end

