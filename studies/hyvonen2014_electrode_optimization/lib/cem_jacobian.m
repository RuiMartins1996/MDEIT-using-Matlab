function J = cem_jacobian(ctx, fwd)
%CEM_JACOBIAN  Adjoint EIT Jacobian J = dV/dsigma, [N*M x n_elem]
%(PLAN_implementation.md sec 3.6).
%
%   J = cem_jacobian(ctx, fwd)     % fwd = cem_fwd_solve(...) output
%
% Row ordering matches V(:) for V = fwd.V ([M x N]): row index = (j-1)*M+i,
% i.e. for pattern j, all M mean-free electrode measurements i=1..M in
% order, then the next pattern. J_{(j,i),k} = -Area_k*grad(y_i).grad(x_j),
% the standard adjoint EIT sensitivity formula (here restricted to the
% "u" (nodal potential) part of the augmented CEM state vector).

n = ctx.n_nodes; M = fwd.M; N = fwd.N;
n_elem = ctx.n_elem;

Xu = fwd.X(1:n,:);   % [n x N]
Yu = fwd.Y(1:n,:);   % [n x M]

J = zeros(N*M, n_elem);
for j = 1:N
    xj = Xu(:,j);
    for i = 1:M
        row = (j-1)*M + i;
        J(row,:) = -local_bilinear_contract(ctx, Yu(:,i), xj)';
    end
end
end
