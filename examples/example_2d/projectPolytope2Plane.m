
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

