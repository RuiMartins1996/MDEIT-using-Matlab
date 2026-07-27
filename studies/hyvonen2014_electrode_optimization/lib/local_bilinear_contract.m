function c = local_bilinear_contract(ctx, u, v)
%LOCAL_BILINEAR_CONTRACT  Per-element bilinear form c(k) = Area_k * grad(u)_k . grad(v)_k
%for two nodal P1 fields u, v (each n_nodes x 1).
%
%   c = local_bilinear_contract(ctx, u, v)    % c is [n_elem x 1]
%
% This is the vectorized building block for both the adjoint EIT
% sensitivity formula (J_{(j,i),k} = -c(k) with u=y_i, v=x_j) and its
% theta-derivative (same formula with (x_j_theta,y_i) or (x_j,y_i_theta)):
% see PLAN_implementation.md sec 3.6 / 4.2 -- both are point evaluations
% of this one routine, so there is only one place a stiffness-contraction
% bug could hide.

elems = ctx.elems;
U = full(u(elems));  % n_elem x 3
V = full(v(elems));

g1 = ctx.elem_grad(:,:,1); g2 = ctx.elem_grad(:,:,2);

gu1 = sum(g1.*U,2); gu2 = sum(g2.*U,2);
gv1 = sum(g1.*V,2); gv2 = sum(g2.*V,2);

c = ctx.elem_area .* (gu1.*gv1 + gu2.*gv2);
end
