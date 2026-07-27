function [x,w] = gauss_legendre_nodes(n)
%GAUSS_LEGENDRE_NODES  Nodes/weights for Gauss-Legendre quadrature on [-1,1].
%
%   [x,w] = gauss_legendre_nodes(n)
%
% Standard Golub-Welsch eigenvalue construction (n small, e.g. 5-40, so
% this is cheap and needs no toolbox).

beta = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
T = diag(beta,1) + diag(beta,-1);
[V,D] = eig(T);
x = diag(D);
[x,idx] = sort(x);
V = V(:,idx);
w = 2*(V(1,:).^2)';
x = x(:);
end
