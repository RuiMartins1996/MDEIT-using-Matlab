function E = electrode_arc_integrals(ctx, theta_minus, theta_plus)
%ELECTRODE_ARC_INTEGRALS  Partial-edge CEM electrode integrals + their
%exact derivatives w.r.t. the arc endpoints (PLAN_implementation.md sec 3.3).
%
%   E = electrode_arc_integrals(ctx, theta_minus, theta_plus)
%
% theta_minus < theta_plus, both real scalars (NOT reduced mod 2pi --
% callers keep a running/unwrapped angle so electrode ordering is
% unambiguous; this function reduces internally).
%
% Computes, for a single electrode occupying the arc [theta_minus,theta_plus]:
%   Bm   = int_{Em} phi_i phi_j ds        [n_nodes x n_nodes] sparse
%   cm   = int_{Em} phi_i ds              [n_nodes x 1] sparse
%   Elen = int_{Em} ds
% and the boundary-term building blocks needed for d(.)/dtheta_minus,
% d(.)/dtheta_plus via the fundamental theorem of calculus:
%   dBm/dtheta_plus  = + speed_plus  * phi_plus  * phi_plus'
%   dBm/dtheta_minus = - speed_minus * phi_minus * phi_minus'
% (and analogously, with phi instead of phi*phi', for cm; and with 1
% instead of phi, for Elen). E.phi_minus, E.phi_plus are sparse n_nodes x 1
% vectors (nonzero only at the 2 nodes of the edge containing that
% endpoint); E.speed_minus, E.speed_plus are |gammadot| there.
%
% Wraparound (arcs crossing theta=0) is handled by testing 3 shifted
% copies (-2pi,0,+2pi) of every boundary edge against the (already
% unwrapped) electrode interval -- see PLAN_implementation.md sec 3.3.

shape = ctx.shape;
n_nodes = ctx.n_nodes;

% Reduce to the base branch [0,2pi) before hunting for overlapping mesh
% edges -- callers (and, in particular, fminunc's internal line search,
% which is not aware of the repulsion barrier and can propose a trial
% theta many revolutions away) may pass theta_minus far outside any
% single 2pi period. Only theta_minus - theta_plus (the width) matters
% for physical placement, so shifting both endpoints by the same integer
% multiple of 2pi is exact, not an approximation.
k = floor(theta_minus/(2*pi));
theta_minus = theta_minus - k*2*pi;
theta_plus  = theta_plus  - k*2*pi;

n1 = ctx.boundary_edges(:,1);
n2 = ctx.boundary_edges(:,2);
t1 = ctx.boundary_edge_theta(:,1);
t2 = ctx.boundary_edge_theta(:,2);
n_edges = numel(n1);

n_quad = 8;
[xg,wg] = gauss_legendre_nodes(n_quad);

Bm_i = []; Bm_j = []; Bm_v = [];
cm_i = []; cm_v = [];
Elen = 0;

shifts = [-2*pi, 0, 2*pi];

for s = 1:3
    sh = shifts(s);
    ta = t1 + sh; tb = t2 + sh;
    lo = max(ta, theta_minus);
    hi = min(tb, theta_plus);
    hits = find(hi > lo);
    if isempty(hits), continue; end

    for r = hits(:)'
        s0 = lo(r); s1 = hi(r);
        edge_len = tb(r)-ta(r);
        % quadrature points in the *shifted* frame (consistent with ta,tb)
        tq = s0 + 0.5*(s1-s0)*(xg+1);
        w_scaled = 0.5*(s1-s0)*wg;

        % local coordinate along the edge (unaffected by the shift)
        ell = (tq - ta(r)) / edge_len;

        % arc-length speed at the *unshifted* physical angle
        [~,~,speed] = boundary_curve(tq - sh, shape);

        phi1 = 1 - ell;   % shape function at n1(r)
        phi2 = ell;       % shape function at n2(r)
        ds = w_scaled .* speed;

        Elen = Elen + sum(ds);

        cm_i(end+1,1) = n1(r); cm_v(end+1,1) = sum(phi1.*ds); %#ok<AGROW>
        cm_i(end+1,1) = n2(r); cm_v(end+1,1) = sum(phi2.*ds); %#ok<AGROW>

        Bm_i(end+1,1) = n1(r); Bm_j(end+1,1) = n1(r); Bm_v(end+1,1) = sum(phi1.*phi1.*ds); %#ok<AGROW>
        Bm_i(end+1,1) = n2(r); Bm_j(end+1,1) = n2(r); Bm_v(end+1,1) = sum(phi2.*phi2.*ds); %#ok<AGROW>
        Bm_i(end+1,1) = n1(r); Bm_j(end+1,1) = n2(r); Bm_v(end+1,1) = sum(phi1.*phi2.*ds); %#ok<AGROW>
        Bm_i(end+1,1) = n2(r); Bm_j(end+1,1) = n1(r); Bm_v(end+1,1) = sum(phi1.*phi2.*ds); %#ok<AGROW>
    end
end

if isempty(Bm_i)
    error('electrode_arc_integrals:noOverlap', ...
        'Electrode arc [%.4f, %.4f] does not overlap any boundary edge -- width too small for this mesh?', ...
        theta_minus, theta_plus);
end

Bm = sparse(Bm_i, Bm_j, Bm_v, n_nodes, n_nodes);
cm = sparse(cm_i, ones(size(cm_i)), cm_v, n_nodes, 1);

% --- boundary terms at the two endpoints, for the FTC derivative -------
[phi_minus, speed_minus] = point_basis(ctx, theta_minus);
[phi_plus,  speed_plus]  = point_basis(ctx, theta_plus);

E = struct();
E.Bm = Bm; E.cm = cm; E.Elen = Elen;
E.phi_minus = phi_minus; E.speed_minus = speed_minus;
E.phi_plus  = phi_plus;  E.speed_plus  = speed_plus;
end

function [phi, speed] = point_basis(ctx, theta)
% Sparse n_nodes x 1 vector of the boundary P1 basis function values at
% angle parameter theta (reduced mod 2pi to locate the containing edge),
% plus the arc-length speed |gammadot(theta)| there.
shape = ctx.shape;
theta_r = mod(theta, 2*pi);

t1 = ctx.boundary_edge_theta(:,1);
t2 = ctx.boundary_edge_theta(:,2);
n1 = ctx.boundary_edges(:,1);
n2 = ctx.boundary_edges(:,2);

shifts = [-2*pi, 0, 2*pi];
found = false;
for s = 1:3
    sh = shifts(s);
    ta = t1 + sh; tb = t2 + sh;
    r = find(theta_r >= ta - 1e-12 & theta_r <= tb + 1e-12, 1);
    if ~isempty(r)
        found = true;
        ell = (theta_r - ta(r)) / (tb(r)-ta(r));
        ell = min(max(ell,0),1);
        phi = sparse([n1(r); n2(r)], [1;1], [1-ell; ell], ctx.n_nodes, 1);
        break
    end
end
if ~found
    error('electrode_arc_integrals:pointNotFound', ...
        'Could not locate boundary edge containing theta=%.6f', theta_r);
end
[~,~,speed] = boundary_curve(theta, shape);
end
