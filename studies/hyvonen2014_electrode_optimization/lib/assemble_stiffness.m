function K = assemble_stiffness(ctx, w)
%ASSEMBLE_STIFFNESS  K(w) = sum_k w_k * K_k, the P1 stiffness matrix with
%per-element weight vector w (n_elem x 1). Used both for the physical
%conductivity (w = sigma) and, in the shape-gradient contraction, for
%arbitrary element-weight vectors (PLAN_implementation.md sec 4.2, where
%w = W_row is an OED-gradient weight, not a conductivity).
%
%   K = assemble_stiffness(ctx, w)

n_elem = ctx.n_elem;
elems = ctx.elems;
local_stiff = ctx.local_stiff;   % [n_elem x 3 x 3]

w = w(:);
ii = zeros(9*n_elem,1); jj = zeros(9*n_elem,1); vv = zeros(9*n_elem,1);
pos = 0;
for a = 1:3
    for b = 1:3
        idx = pos*n_elem+1 : (pos+1)*n_elem;
        ii(idx) = elems(:,a);
        jj(idx) = elems(:,b);
        vv(idx) = w .* squeeze(local_stiff(:,a,b));
        pos = pos+1;
    end
end
K = sparse(ii,jj,vv,ctx.n_nodes,ctx.n_nodes);
end
