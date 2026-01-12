function mesh_quality_metrics = compute_mesh_quality_metrics(fmdl)
    
    mesh_quality_metrics = struct();
    
    aspect_ratio_per_element = compute_aspect_ratio(fmdl);
    radius_ratio_per_element = compute_radius_ratio(fmdl);
    mean_ratio_per_element = compute_mean_ratio(fmdl);
    shape_metric_per_element = compute_jacobian_condition_number(fmdl);

    % % DEBUG:
    % regular.elems = [1,2,3,4];
    % regular.nodes = [0 0 0;...
    %     1,0,0;...
    %     1/2,sqrt(3)/2,0;...
    %     1/2,sqrt(3)/6,sqrt(6)/3];
    % 
    % irregular.elems = [1,2,3,4];
    % irregular.nodes = [0 0 0;...
    %     1,0,0;...
    %     0,1,0;...
    %     0,0,1];

    mesh_quality_metrics.aspect_ratio_per_elemet = aspect_ratio_per_element;
    mesh_quality_metrics.radius_ratio_per_element = radius_ratio_per_element;
    mesh_quality_metrics.mean_ratio_per_element = mean_ratio_per_element;
    mesh_quality_metrics.shape_metric_per_element = shape_metric_per_element;
end


function aspect_ratio_per_element = compute_aspect_ratio(fmdl)
nodes = fmdl.nodes;
elems = fmdl.elems;

% Extract node coordinates for each element
X1 = nodes(elems(:,1), :);
X2 = nodes(elems(:,2), :);
X3 = nodes(elems(:,3), :);
X4 = nodes(elems(:,4), :);

% Compute all 6 edge lengths
l12 = sqrt(sum((X1 - X2).^2, 2));
l13 = sqrt(sum((X1 - X3).^2, 2));
l14 = sqrt(sum((X1 - X4).^2, 2));
l23 = sqrt(sum((X2 - X3).^2, 2));
l24 = sqrt(sum((X2 - X4).^2, 2));
l34 = sqrt(sum((X3 - X4).^2, 2));

% Combine into a single matrix of edge lenghts
edge_lengths = [l12, l13, l14, l23, l24, l34];

% Compute mesh aspect ratio
aspect_ratio_per_element = max(edge_lengths, [], 2) ./ min(edge_lengths, [], 2);
end

function Q_RR = compute_radius_ratio(fmdl)
% COMPUTE_RADIUS_RATIO  Compute the ratio between the in-radius and circumradius 
% mesh quality metric for tetrahedra

    nodes = fmdl.nodes;
    elems = fmdl.elems;
    Ne = size(elems, 1);

    Q_RR = zeros(Ne, 1);

    for k = 1:Ne
        % Get node coordinates
        X = nodes(elems(k,:), :);
        x1 = X(1,:)'; x2 = X(2,:)'; x3 = X(3,:)'; x4 = X(4,:)';

        % --- Volume ---
        V = abs(dot(x2 - x1, cross(x3 - x1, x4 - x1))) / 6;

        % --- Face areas ---
        A1 = 0.5 * norm(cross(x2 - x3, x4 - x3)); % opposite x1
        A2 = 0.5 * norm(cross(x1 - x3, x4 - x3)); % opposite x2
        A3 = 0.5 * norm(cross(x1 - x2, x4 - x2)); % opposite x3
        A4 = 0.5 * norm(cross(x1 - x2, x3 - x2)); % opposite x4

        % --- Inradius --- (https://en.wikipedia.org/wiki/Tetrahedron)
        r_in = 3 * V / (A1 + A2 + A3 + A4);

        % --- Circumradius ---
        % Construct matrix of relative edge vectors
        a = x2 - x1;
        b = x3 - x1;
        c = x4 - x1;
        M = [a b c];

        % Right-hand side for circumcenter formula
        rhs = [dot(x2,x2)-dot(x1,x1); dot(x3,x3)-dot(x1,x1); dot(x4,x4)-dot(x1,x1)] / 2;

        % Solve for circumcenter in local coordinates
        C = M' \ rhs;
        R_circ = norm(C-x1);

        % --- Radius ratio ---
        Q_RR(k) = 3 * r_in / R_circ;
    end
end

function  shape_metric_per_element = compute_jacobian_condition_number(fmdl)
% COMPUTE_JACOBIAN_CONDITION  Compute condition number of Jacobian for each tetrahedron
% and the shape metric

% In Tetrahedral Element Shape Optimization via the Jacobian Determinant and Condition Number

nodes = fmdl.nodes;
elems = fmdl.elems;
Ne = size(elems,1);

shape_metric_per_element = zeros(Ne,1);

W = [1,0,0;...
     1/2,sqrt(3)/2,0;...
     1/2,sqrt(3)/6,sqrt(6)/3]'; % Jacobian matrix of the logical tetrahedron

% Compute  weighted condition number of the matrix An
for k = 1:Ne
    % Node coordinates
    X = nodes(elems(k,:), :);
    x1 = X(1,:)'; x2 = X(2,:)'; x3 = X(3,:)'; x4 = X(4,:)';

    % Construct Jacobian
    J = [x2 - x1, x3 - x1, x4 - x1];

    % Check for degenerate element
    if abs(det(J)) < 1e-12
        Q_CN(k) = Inf; % degenerate element
        continue;
    end
    
    % Condition number
    kappa = norm(J/W,'fro')*norm(W/J,'fro');

    shape_metric_per_element(k) = 3/kappa;
end

end

function Q_MR = compute_mean_ratio(fmdl)
% COMPUTE_MEAN_RATIO  Compute mean ratio (volume-to-edge-length ratio) for tetrahedra

%Joe–Liu quality measure, eq.(4) in:
% A method for reconstructing tomographic images
% of evoked neural activity with electrical impedance
% tomography using intracranial planar arrays

nodes = fmdl.nodes;
elems = fmdl.elems;
Ne = size(elems,1);

Q_MR = zeros(Ne,1);

for k = 1:Ne
    % Node coordinates
    X = nodes(elems(k,:), :);
    x1 = X(1,:)'; x2 = X(2,:)'; x3 = X(3,:)'; x4 = X(4,:)';

    % --- Volume ---
    V = abs(dot(x2 - x1, cross(x3 - x1, x4 - x1))) / 6;

    % --- Edge lengths ---
    l12 = norm(x1-x2); l13 = norm(x1-x3); l14 = norm(x1-x4);
    l23 = norm(x2-x3); l24 = norm(x2-x4); l34 = norm(x3-x4);

    L2_sum = l12^2 + l13^2 + l14^2 + l23^2 + l24^2 + l34^2;

    % --- Mean ratio ---
    Q_MR(k) = (12 * (3*V)^(2/3)) / L2_sum;
end
end
