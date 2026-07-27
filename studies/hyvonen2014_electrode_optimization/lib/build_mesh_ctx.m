function ctx = build_mesh_ctx(shape, maxsz)
%BUILD_MESH_CTX  Precompute a fixed P1 triangular mesh + geometry for a
%star-shaped 2D domain, once, for reuse across every electrode
%configuration (Deviation 2 in PLAN_implementation.md: we never remesh
%inside the optimization loop).
%
%   ctx = build_mesh_ctx(shape, maxsz)
%
% shape.type = 'disk', shape.radius = R. maxsz is Netgen's target element
% size (passed straight to ng_mk_cyl_models).
%
% ctx fields:
%   nodes        [n_nodes x 2]
%   elems        [n_elem  x 3]  node indices (P1 triangles)
%   n_nodes, n_elem
%   elem_area    [n_elem x 1]
%   elem_grad    [n_elem x 3 x 2]  gradient of each of the 3 local P1
%                shape functions (constant per element)
%   local_stiff  [n_elem x 3 x 3]  element stiffness contribution with
%                unit conductivity (Area * grad_a . grad_b)
%   boundary_edges      [n_edges x 2]  node indices, oriented so that
%                        going node1->node2 increases the angle parameter
%   boundary_edge_theta  [n_edges x 2] unwrapped angle parameter of the
%                        two edge endpoints (theta2 > theta1, both in a
%                        single branch, NOT reduced mod 2pi)
%   node_theta   [n_nodes x 1]  atan2(y,x) reduced to [0,2pi) (only
%                meaningful/used for boundary nodes)
%   shape        (echoed back)

fmdl = ng_mk_cyl_models([0, shape.radius, maxsz], [0], []);

nodes = fmdl.nodes;
elems = fmdl.elems;

n_nodes = size(nodes,1);
n_elem  = size(elems,1);

% --- Per-element P1 geometry -------------------------------------------------
elem_grad   = zeros(n_elem,3,2);
elem_area   = zeros(n_elem,1);
local_stiff = zeros(n_elem,3,3);

for k = 1:n_elem
    idx = elems(k,:);
    P = nodes(idx,:);   % 3x2
    % signed area x2 via cross product
    v1 = P(2,:)-P(1,:); v2 = P(3,:)-P(1,:);
    twoA = v1(1)*v2(2) - v1(2)*v2(1);
    A = abs(twoA)/2;
    elem_area(k) = A;

    % gradients of P1 basis functions: grad(phi_i) = perp(edge opposite i)/(2A)
    % standard formula using rotated edge vectors
    e = [P(2,:)-P(3,:); P(3,:)-P(1,:); P(1,:)-P(2,:)]; % edges opposite node 1,2,3
    % rotate by 90deg: (x,y) -> (y,-x), then sign so that grad points "inward" consistently
    rot = [e(:,2), -e(:,1)];
    grad = rot / twoA;   % using signed 2A preserves correct orientation regardless of winding
    elem_grad(k,:,:) = reshape(grad,1,3,2);

    local_stiff(k,:,:) = reshape(A * (grad*grad'),1,3,3);
end

% --- Boundary edges -----------------------------------------------------
% An edge belongs to the boundary iff it appears in exactly one triangle.
E = [elems(:,[1 2]); elems(:,[2 3]); elems(:,[3 1])];
Es = sort(E,2);
[Eu,~,ic] = unique(Es,'rows');
counts = accumarray(ic,1);
is_boundary = counts(ic) == 1;
boundary_pairs_raw = E(is_boundary,:);   % as originally listed (some orientation)

theta_all = mod(atan2(nodes(:,2),nodes(:,1)), 2*pi);

n1 = boundary_pairs_raw(:,1);
n2 = boundary_pairs_raw(:,2);
t1 = theta_all(n1);
t2 = theta_all(n2);

% signed shortest angular difference n1->n2, wrapped to (-pi,pi]
dtheta = mod(t2 - t1 + pi, 2*pi) - pi;

% orient so that dtheta > 0 (n1 -> n2 goes counter-clockwise / increasing theta)
swap = dtheta < 0;
tmp = n1; n1(swap) = n2(swap); n2(swap) = tmp(swap);
tmp = t1; t1(swap) = t2(swap); t2(swap) = tmp(swap);
dtheta(swap) = -dtheta(swap);

theta1_unwrapped = t1;
theta2_unwrapped = t1 + dtheta;   % > theta1_unwrapped, small positive gap

ctx = struct();
ctx.nodes = nodes;
ctx.elems = elems;
ctx.n_nodes = n_nodes;
ctx.n_elem  = n_elem;
ctx.elem_area = elem_area;
ctx.elem_grad = elem_grad;
ctx.local_stiff = local_stiff;
ctx.boundary_edges = [n1, n2];
ctx.boundary_edge_theta = [theta1_unwrapped, theta2_unwrapped];
ctx.node_theta = theta_all;
ctx.shape = shape;

fprintf('Mesh (2D disk, Hyvonen 2014 CEM): %i elements, %i nodes, %i boundary edges\n',...
    n_elem, n_nodes, size(ctx.boundary_edges,1));
end
