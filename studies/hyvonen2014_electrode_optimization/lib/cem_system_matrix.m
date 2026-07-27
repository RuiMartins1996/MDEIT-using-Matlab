function [A_red, M_sel] = cem_system_matrix(ctx, sigma, elec)
%CEM_SYSTEM_MATRIX  Grounded CEM system matrix (PLAN_implementation.md sec 3.3-3.4).
%
%   [A_red, M_sel] = cem_system_matrix(ctx, sigma, elec)
%
% elec = assemble_electrodes(...) output, M = numel(elec) electrodes.
% Unknowns x_red = [u (n_nodes); U_1; ...; U_{M-1}], with electrode M
% GROUNDED (U_M := 0, fixed). This is a valid gauge choice for the CEM
% (which is only defined up to an additive constant on (u,U)): deleting
% the last row+column of the (n+M)x(n+M) system is standard "ground a
% reference electrode" practice, and is exact here because the dropped
% row is linearly dependent on the others whenever the injected currents
% sum to zero (which our patterns always do).
%
% M_sel: [M x (n+M-1)] so that V = M_sel*x_red gives the MEAN-FREE vector
% of all M electrode potentials (U_1,...,U_{M-1},0) for this x_red -- this
% is the representative of U(sigma) in C^M/C that all downstream
% quantities (data, Jacobian, adjoint) use consistently.

n = ctx.n_nodes;
M = numel(elec);

K = assemble_stiffness(ctx, sigma);
KB = K;
for m = 1:M
    KB = KB + elec(m).Bm/elec(m).z;
end

Cblock = sparse(n, M-1);
Ddiag  = zeros(M-1,1);
for m = 1:M-1
    Cblock(:,m) = -elec(m).cm/elec(m).z;
    Ddiag(m)    = elec(m).Elen/elec(m).z;
end
D = spdiags(Ddiag,0,M-1,M-1);

A_red = [KB, Cblock; Cblock.', D];
A_red = (A_red+A_red.')/2;   % force exact symmetry (roundoff)

P0 = eye(M) - ones(M,M)/M;
M_sel = sparse(M, n+M-1);
M_sel(:, n+1:n+M-1) = P0(:,1:M-1);
end
