function [xk,k,allxk] = gaussNewtonV2(res,jac,x0,tol,maxIterations)

% The damping parameter is going to be adaptive
lambda = 0.1;      
lambda_inc = 10;    % factor to increase if step rejected
lambda_dec = 0.1;   % factor to decrease if step accepted

lambda_min = 1e-12; % safeguard
lambda_max = 1e+6;  % safeguard

% Initial condition 
xk = x0;

allxk = nan(length(x0),maxIterations+1);
allxk(:,1) = x0;

rk = res(x0);
normRk = norm(r0,2);

% The objective function \phi(\sigma) = 0.5*\sum_i r_i^2(\sigma)
phi_old = 0.5*norm(rk,2)^2;

k= 0;
while true
    
    % Print diagonostics
    fprintf('gn | iteration %i | (r = %.5f)\n',k,normRk)

    % Compute resiudal
    rk = res(xk);
    normRk = norm(rk,2);

    % Stopping criterion
    if normRk<tol
        fprintf('Terminated because ||rk||< tol \n');
        break
    end
    
    % Stability check
    if normR0<normRk
        error('Bad step. Diverging!');
    end
    
    % Compute Jacobian
    Jk = jac(xk);
    
    % Construct NOSER prior
    R = sparse(1:length(x0),1:length(x0),sum(Jk.^2,1),length(x0),length(x0));
    
    % lhs matrix
    A = Jk'*Jk+lambda*R;
    b = Jk'*-rk;
    
    %Simple Jacobi preconditioner for the pcg method 
    L = sqrt(diag(diag(A)));
    U = L;
    
    % Compute search direction by solving normal equations
    pk = pcg(A,b,1e-5,numel(x0),L,U);

    % Trial step
    x_trial = xk + pk;
    rk_trial = res(x_trial);
    phi_new = 0.5*norm(rk_trial,2)^2;

    % Approximate predicted decrease
    pred  = 0.5 * pk' * (lambda*pk - Jk'*rk);














    % Update solution
    xk = xk + pk;

    % Update number of iterations
    k = k + 1;
    
    % Store vector
    allxk(:,k+1) = xk;

    if k>maxIterations
        normRk = norm(res(xk),2);
        fprintf('gn | iteration %i | (r = %.5f)\n',k,normRk)
        fprintf('gn exceeded max # of iterations: %i\n',maxIterations);
        break;
    end
end

end

