function [U,S,V] = rsvd(A,k,tol)
% Implementation of randomized svd decomposition: https://gregorygundersen.com/blog/2019/01/17/randomized-svd/

% INPUTS:


% tol: tolerance for adaptive range finder

% k : target rank

% l : sampling parameter indicating the number of Gaussian random vectors
% to draw for Omega

m = size(A,1);
n = size(A,2);

use_adaptive_range_finder = false;

if nargin<3
    tol = 1e-6;
end

if nargin>2
    use_adaptive_range_finder = true;
end

% For Gaussian test matrices, it is adequate to choose the oversampling parameter to be a small constant, 
% such as p=5 or p=10. There is rarely any advantage to select p>k.
p = 5;

if use_adaptive_range_finder
    % Algorithm 3: Algorithm 3: Iterative construction of Q
    maxit = min(size(A));
    Q = adaptive_range_finder(A, tol, maxit);
else
    % % Algorithm 2: Randomized range finder
    Omega = randn(n, k+p);
    Y = A * Omega; % [m x l]
    [Q,~] = qr(Y,0); %QR factorization of Y to generate an orthonormal matrix Q
end


% Algorithm 1. Randomized SVD
B = Q' * A;
[Uh,S,V] = svd(B,'econ');

U = Q * Uh;

end


function Q = adaptive_range_finder(A, tol, maxit)
% Adaptive randomized range finder
% A     : matrix or function handle for A*x
% tol   : relative tolerance
% maxit : safety cap on iterations

[m, ~] = size(A);
Q = zeros(m, maxit);  % preallocate
col = 0;

Omega = randn(size(A,2),maxit);
for i = 1:maxit

    % 1. Random probe
    omega = Omega(:,i);

    % 2. Apply A
    y = A * omega;

    % 3. Orthogonalize
    if col > 0
        y = y - Q(:,1:col)*(Q(:,1:col)'*y);  % orthogonalize
        y = y - Q(:,1:col)*(Q(:,1:col)'*y);  % optional double pass
    end

    % 4. Normalize
    ny = norm(y);
    if ny < 1e-14
        continue
    end
    q = y / ny;

    % 5. Expand basis
    col = col + 1;

    Q(:,col) = q;

    % 6. Check residual norm (cheap estimator)
    if mod(i,100) == 0
        residual = estimate_residual(A, Q(:,1:col));
        fprintf('Iter %i/%i| r = %2.2g, tol = %2.2g\n',i,maxit,residual,tol);
        if residual < tol
            break
        end
    end
end

Q = Q(:,1:col);  % trim at the end

end


function err = estimate_residual(A, Q)
nr = 3;  % number of random probes
err = 0;

G = randn(size(A,2),nr);

for j = 1:nr
    r = A*G(:,j) - Q*(Q'*(A*G(:,j)));
    err = max(err, norm(r));
end
end