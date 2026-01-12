function [xk,num_of_iterations,pcg_iterations,residual,exit_flag] = ...
    gn_solve(res,jac,x0,tol,max_iterations,rt_r,lambda,M,verbose)

if nargin < 9
    verbose = false;
end

if nargin < 8
    M = [];
end

if not(isa(rt_r, 'function_handle'))
    error('RtR should be a function handle')
end

min_num_of_steps = ceil(max_iterations/10);

% Initial state
xk = x0;

% Compute residual and jacobian
rk = res(xk);
jk = jac(xk);

% Compute objective and gradient at xk
fk1 = 0.5*(rk'*rk);  %avoid calling function func to not recompute res(xk)
gk = jk'*rk;        %avoid calling function grad_func to not recompute res(xk) and jac(xk)


num_of_iterations = 0;
pcg_iterations = [];

% Exit imediatly if convergence criterion is satisfied
if norm(gk,'inf')<tol
    residual = rk;
    exit_flag = 3;
    fprintf('Critical point at ||g|| =  %.2g\n',norm(gk));
    fprintf('GN time: %.2f\n', toc);
    return;
end

if verbose
report_diagonostic_info(0,norm(rk,2),norm(jac(x0)'*rk))
end

tic; 
while true

    %rk,jk,fk and gk have already been computed
    
    % Matrix-free operators ( assemble A and b only once outside pcg)
    Afun  = @(x) jk' * (jk*x) + lambda * rt_r(x,jk);
    b = jk' * (-rk);

    % Solve Afun * p = b (use pcg with optional preconditioner M)
    if isempty(M)
        [pk,flag,relres,iter] = pcg(Afun, b, tol, numel(x0));
        pcg_iterations(num_of_iterations+1) = iter;

    else
        [pk,flag,relres,iter] = pcg(Afun, b, tol, numel(x0), M);
        pcg_iterations(num_of_iterations+1) = iter;
    end
    
    % Accept the step
    xk = xk + pk;

    % Update counter
    num_of_iterations = num_of_iterations+1;

    % Report diagonostic info
    if verbose
        report_diagonostic_info(num_of_iterations,norm(rk,2),norm(-b))
    end
    
    rk1 = rk;
    % Update residual
    rk = res(xk);
    
    if norm(rk,'inf')>norm(rk1,'inf')
        residual = rk1;
        xk = xk-pk;

        exit_flag = 6;
        fprintf('GN time: %.2f\n', toc);
        fprintf('Residual increase');
        return;
    end

    % Stopping Criterion (1) : maximum number of iterations
    if num_of_iterations >= max_iterations
        residual = rk;
        exit_flag = 1;
        fprintf('GN time: %.2f\n', toc);
        fprintf('GN exceeded max # of iterations: %i\n', max_iterations);
        return;
    end

    % Stopping criterion (4): Minimum number of steps, but no change in solution
    if num_of_iterations>min_num_of_steps && norm(pk)<=tol*(norm(xk)+tol)
        residual = rk;
        exit_flag = 4;
        fprintf('No change in solution ||pk|| = %.2g\n',norm(pk));
        fprintf('GN time: %.2f\n', toc);
        return;
    end
    
    % Update Jacobian
    jk = jac(xk);

    % Compute objective and gradient at xk
    fk = 0.5*(rk'*rk);  %avoid calling function func to not recompute res(xk)
    gk = jk'*rk;        %avoid calling function grad_func to not recompute res(xk) and jac(xk)

    
    % Stopping criterion (3): Gradient of cost function is small
    if norm(gk,'inf')<tol
        residual = rk;
        exit_flag = 3;
        fprintf('Critical point at ||g|| =  %.2g\n',norm(gk));
        fprintf('GN time: %.2f\n', toc);
        return;
    end

    % Stopping criterion (5): Cost function no longer decreasing
    if abs(fk1-fk)<=tol && num_of_iterations>min_num_of_steps
        residual = rk;
        exit_flag = 5;
        fprintf('No change in cost function\n');
        fprintf('GN time: %.2f\n', toc);
        return
    end
    
    % Update cost function value at last iteration
    fk1 = fk;


end

end



function report_diagonostic_info(iteration,residual_norm,grad_norm)

fprintf('GN | it %i | r=%.5g| g=%.5g\n',...
    iteration,residual_norm,grad_norm);
end


function out = func(x,res)
    r = res(x);
    out = 0.5*(r'*r);
end

function out = grad_func(x,res,jac)
    r = res(x);
    j = jac(x);
    out = j'*r;
end