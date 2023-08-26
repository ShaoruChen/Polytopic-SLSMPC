% Plot the region of attraction (feasible domain) of each robust MPC
% method with given uncertainty parameters. 

clear
close all

%% eps_A comparison
data = load('data/example_2d_diags_list_eps_A_hor_3.mat');
diags_list_eps_A = data.diags_list;
eps_A_list = 0.05:0.05:0.45;

N_eps = length(eps_A_list);
method_cell = {'Tube_MPC', 'Tube_MPC_Homothetic', 'Tube_MPC_Nominal', 'Tube_MPC_Flexible',  ...
               'SLS_MPC', 'Lumped_Dist_MPC', 'Constr_Tightening_MPC'};
num_method = 7;

ii= 8;

diags_record = diags_list_eps_A{ii};
N = length(diags_record);

% find the best coverage for each tube-based MPC method
best_coverage = zeros(1, length(method_cell));
feasible_x0 = cell(1, length(method_cell));

for jj = 1:N
    diags = diags_record{jj};
    method = diags.method;
    ind = find(strcmp(method_cell, method));
    
    coverage = diags.feasible_rate;
    if coverage >= best_coverage(ind)
        best_coverage(ind) = coverage;
        feasible_x0{ind} = diags.feasible_set;
    end
end


data = load('data/MPC_example_data.mat');
MPC_data = data.MPC_data;

%% plot the feasible regions

figure;
hold on
P_reduced = projectPolytope2Plane(MPC_data.terminal_constr);
ps = polyshape(P_reduced.V(:,1), P_reduced.V(:,2));
pg = plot(ps);
pg.FaceColor = 'none';
pg.LineWidth = 2.0;

x0_set = diags_record{1}.x0_set;
scatter(x0_set(:,1), x0_set(:,2), 'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 0.7*[1 1 1]);

method_cell = {'Tube-A', 'Tube-B', 'Tube-C', 'Tube-D',  ...
               'SLS-MPC', 'Lumped-Disturbance', 'Offline-Tightening'};

colors = {"#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"};

h = {};
count = 1;
for ii = 1:length(method_cell)
    feasible_x = feasible_x0{ii};
    if ~isempty(feasible_x)
        feasible_set = Polyhedron(feasible_x);
        feasible_set.minHRep();
        feasible_set.minVRep();
        P_reduced = projectPolytope2Plane(feasible_set);
        ps = polyshape(P_reduced.V(:,1), P_reduced.V(:,2));
        h{count} = plot(ps);
        h{count}.FaceColor = 'none';
        h{count}.EdgeColor = colors{ii};
        h{count}.LineWidth = 2.0;
        h{count}.DisplayName = method_cell{ii};

        if ii <= 4
           h{count}.LineStyle = '-.';
        end
        count = count + 1;
    end
end

xlim([-3.3 3.3]);
ylim([-8.5 8.5]);

ax = gca; % get handle to the current axes
ax.FontSize = 18;
% lg = legend(ax, h{:}, 'Orientation', 'Horizontal');
legend([h{:}], 'Location', 'bestoutside');


xlabel('$x_1$', 'Interpreter', 'Latex', 'FontSize', 24);
ylabel('$x_2$', 'Interpreter', 'Latex', 'FontSize', 24);

