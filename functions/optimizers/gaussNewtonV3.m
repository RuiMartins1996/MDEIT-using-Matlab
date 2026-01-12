function [xk,numOfIterations,allxk] = gaussNewtonV3(res,jac,x0,tol,maxIterations,Wsqrt,RtR,lambda)

if not(isa(RtR, 'function_handle'))
    error('RtR should be a function handle')
end

if isempty(Wsqrt)
    Wsqrt = 1;
end

if nargin < 8 || isempty(lambda) 
    lambda = 1;
end

xk = x0;

allxk = nan(length(x0),maxIterations+1);
allxk(:,1) = x0;

r0  = res(x0);
r0 = Wsqrt*r0;
normR0 = norm(r0,2);

fprintf('LM | it %i | r=%.5f \n',0,normR0);

numOfIterations = 0;

% Compute Jacobian
Jk = jac(xk);
Jk = Wsqrt*Jk;

% This should be updated later. The rank does not need to be computed like
% this. We know the rank as a function of the number of measurements
rankJk = rank(Jk); % Rank of the Jacobian matrix

% Discrepancy principle (don't need delta if using prewhitened residual)
tau = 1.5;

rk = r0;

while true
    
    % Since J^T*J is symmetric, then in the eigenvalue decomposition V*D*V^T,
    % it's true that V^T*V = 1

    % The singular value decomposition of J is U S V^T, where V is the same
    % as above. 

    [U,S,V] = svds(Jk,rankJk,'largest');
    
    %this comes from working out the GN normal equations with the SVD for J 
    pk = -V*((U'*rk)./diag(S)); 
    
    % Use regularization to avoid amplifying noise
    % pk = linsolve(V*S.^2*V'+lambda*RtR(Jk),-V*S*U'*rk);

    % Update solution
    xk = xk + pk;
   
    % Compute resiudal
    rk = res(xk);
    rk = Wsqrt*rk;
    
    %Stability check
    if norm(rk,2)>=normR0
        error('Bad step!')
    end
    
    % 
    if -pk'*(Jk'*Jk)*pk>0
        error('pk is not descent direction');
    end

    % Compute Jacobian
    Jk = jac(xk);
    Jk = Wsqrt*Jk;
    
    % Report diagonostic info
    numOfIterations = numOfIterations+1;
    fprintf('LM | it %i | r=%.5f \n',numOfIterations,norm(rk,2));

    allxk(:,numOfIterations+1) = xk;

    % Stopping criterion (discrepancy principle applied to the
    % pre-whitened residual)
    if norm(Wsqrt*rk,2) <= tau*sqrt(length(rk)) && numOfIterations>=maxIterations
        fprintf('Converged at iteration %i with residual %.5e\n', ...
            numOfIterations, norm(rk,2));
        break;
    elseif numOfIterations>=maxIterations
        fprintf('LM exceeded max # of iterations: %i\n',maxIterations);
        break
    end
end

end

