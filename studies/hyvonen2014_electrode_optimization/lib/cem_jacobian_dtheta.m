function dJ = cem_jacobian_dtheta(ctx, fwd, elec)
%CEM_JACOBIAN_DTHETA  dJ/dtheta_minus(m) for every design variable m
%(PLAN_implementation.md sec 4.2 -- simplified relative to the plan's
%original per-(j,i) contraction scheme: since M (design variables) is
%small, we form the full [N*M x n_elem] dJ/dtheta_m matrix for each m
%directly, exactly mirroring how example_anomaly_circle_2d.m's grad_z
%forms dJ/dp per sensor-position component. This needs only M+N extra
%solves per design variable (state derivatives x_j_theta, y_i_theta),
%reusing the ONE factorization already computed in fwd.dec -- not the
%O(N*M) solves the plan's original write-up estimated.
%
%   dJ = cem_jacobian_dtheta(ctx, fwd, elec)   % dJ{m} is [N*M x n_elem]
%
% Derivation: x^(j) solves A(theta) x^(j) = b^(j) with b^(j) independent
% of theta (currents are assigned to electrode INDICES, not positions),
% so dx^(j)/dtheta_m = -A^-1 (dA_m x^(j)). Likewise y_i = A^-1 m_i with m_i
% (the measurement selector row) independent of theta, so
% dy_i/dtheta_m = -A^-1 (dA_m y_i). Then, since
% J_{(j,i),k} = -Area_k*grad(y_i)_k . grad(x_j)_k  (cem_jacobian.m),
% product rule gives
% dJ_{(j,i),k}/dtheta_m = -Area_k*[ grad(dy_i/dtheta_m)_k.grad(x_j)_k
%                                  + grad(y_i)_k.grad(dx_j/dtheta_m)_k ].

n = ctx.n_nodes; M = fwd.M; N = fwd.N; n_elem = ctx.n_elem;
dec = fwd.dec;

Xu = fwd.X(1:n,:);
Yu = fwd.Y(1:n,:);

dJ = cell(1,M);
for m = 1:M
    dAm = cem_dA_dthetaminus(ctx, elec, m);

    Xtheta = dec \ full(-(dAm*fwd.X));   % [(n+M-1) x N]
    Ytheta = dec \ full(-(dAm*fwd.Y));   % [(n+M-1) x M]

    Xtheta_u = Xtheta(1:n,:);
    Ytheta_u = Ytheta(1:n,:);

    dJm = zeros(N*M, n_elem);
    for j = 1:N
        xj = Xu(:,j); xjt = Xtheta_u(:,j);
        for i = 1:M
            row = (j-1)*M + i;
            yi = Yu(:,i); yit = Ytheta_u(:,i);
            dJm(row,:) = -( local_bilinear_contract(ctx, yit, xj) ...
                          + local_bilinear_contract(ctx, yi, xjt) )';
        end
    end
    dJ{m} = dJm;
end
end
