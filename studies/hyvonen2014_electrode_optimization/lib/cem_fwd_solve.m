function out = cem_fwd_solve(ctx, sigma, elec, I)
%CEM_FWD_SOLVE  Forward + adjoint solves for the grounded CEM system.
%
%   out = cem_fwd_solve(ctx, sigma, elec, I)
%
% I: [N x M] current patterns (rows = build_current_patterns(M) output).
%
% out.A_red, out.M_sel : see cem_system_matrix.m
% out.dec       : decomposition(A_red) — the ONE factorization, reused for
%                 every forward solve, adjoint solve, and (by the caller)
%                 every shape-derivative solve.
% out.X         : [(n+M-1) x N]  forward states, X(:,j) = A_red \ b_red^(j)
% out.Y         : [(n+M-1) x M]  adjoint states, Y(:,i) = A_red \ M_sel(i,:)'
% out.V         : [M x N]  mean-free electrode-potential data, V = M_sel*X

[A_red, M_sel] = cem_system_matrix(ctx, sigma, elec);
n = ctx.n_nodes; M = numel(elec); N = size(I,1);

dec = decomposition(A_red,'ldl');

B = zeros(n+M-1, N);
for j = 1:N
    B(n+1:n+M-1, j) = I(j,1:M-1)';
end
X = dec \ B;

Y = dec \ full(M_sel)';

V = M_sel*X;

out = struct('A_red',A_red,'M_sel',M_sel,'dec',dec,'X',X,'Y',Y,'V',V, ...
    'n_nodes',n,'M',M,'N',N);
end
