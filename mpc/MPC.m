classdef MPC < handle
    %MPC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = public)
        system; % uncertain LTI system: UncertainLTISystem
        horizon;  % horizon used for MPC
        Q; R; Q_T; % quadratic stage cost
        state_constr; % state constraints: Polyhedron
        input_constr; % input constraints: Polyhedron
        terminal_constr; % terminal constraints: Polyhedron
        nx; nu; % state and input dimensions
        x0; % initial condition
    end
    
    methods (Access = public)
        function obj = MPC(params)
            %MPC Construct an instance of this class
            %   Detailed explanation goes here
            obj.system = params.system;
            obj.state_constr = params.state_constr;
            obj.input_constr = params.input_constr;
            obj.terminal_constr = params.terminal_constr;
            obj.horizon = params.horizon;
            obj.Q = params.Q; 
            obj.Q_T = params.Q_T;
            obj.R = params.R;
            obj.nx = obj.system.nx; obj.nu = obj.system.nu;
            if isfield(params, 'x0')
               obj.x0 = params.x0;
            else
               obj.x0 = zeros(obj.nx, 1);
            end
        end
    end
    
    methods (Abstract)
        % solve the finite-horizon robust optimal control problem
        solve(obj) 
    end
    
end

