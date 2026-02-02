function [U,S,V] = rsvd(A,k)
% Implementation of randomized svd decomposition: https://gregorygundersen.com/blog/2019/01/17/randomized-svd/

% INPUTS:

% k : target rank

% l : sampling parameter indicating the number of Gaussian random vectors
% to draw for Omega

m = size(A,1);
n = size(A,2);

if nargin<2
    use_adaptive_range_finder = false;
    k = min(m,n);       % target rank
else
    use_adaptive_range_finder = true;
end

% For Gaussian test matrices, it is adequate to choose the oversampling parameter to be a small constant, 
% such as p=5 or p=10. There is rarely any advantage to select p>k.
p = 10;

if use_adaptive_range_finder
    % Algorithm 3: Algorithm 3: Iterative construction of Q
    Q = adaptive_range_finder(A, 1e-5, 20);
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

[m,~] = size(A);
Q = zeros(m,0);

for i = 1:maxit
    % 1. Random probe
    omega = randn(size(A,2),1);

    % 2. Apply A
    y = A * omega;

    % 3. Orthogonalize
    if ~isempty(Q)
        y = y - Q*(Q'*y);
    end

    % 4. Normalize
    ny = norm(y);
    if ny < 1e-14
        continue
    end
    q = y / ny;

    % 5. Expand basis
    Q = [Q, q];

    % 6. Check residual norm (cheap estimator)
    if estimate_residual(A, Q) < tol
        break
    end
end
end


function err = estimate_residual(A, Q)
nr = 3;  % number of random probes
err = 0;

for j = 1:nr
    g = randn(size(A,2),1);
    r = A*g - Q*(Q'*(A*g));
    err = max(err, norm(r));
end
end