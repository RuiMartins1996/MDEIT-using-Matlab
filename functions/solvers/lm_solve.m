function [xk,num_of_iterations,pcg_iterations,residual,exit_flag] = ...
    lm_solve(res,jac,x0,tol,max_iterations,rt_r,lambda0,M,verbose)
% Implementation of the levenberg-marquadt method described in 
% METHODS FOR NON-LINEAR LEAST SQUARES PROBLEMS, K.Madsen,H.B.Nielsen,O.Tingleff

nu = 2;
% max_num_of_failed_steps = max(1,floor(max_iterations/4));
max_num_of_failed_steps = 20;
min_num_of_steps = max(1,floor(max_iterations/4));

f = @(x) func(x,res);
grad = @(x) grad_func(x,res,jac);

%  Initialize parameter structure for StepSizeSW function call.
params = struct('ftol', 1e-4, 'gtol', 0.9, 'xtol', 1e-6, 'stpmin', 0, ...
    'stpmax', 1e20, 'maxfev', 10000);

if nargin < 9
    verbose = false;
end

if nargin < 8
    M = [];
end

% rt_r is an implementation of direct matrix-vector product (matrix is not
% assembled)
if ~isa(rt_r, 'function_handle')
    error('rt_r must be a function handle');
end

% Initial lambda
lambda = lambda0;

% Initial state
xk  = x0;

% Jacobian and residual at xk
jk = jac(xk);
rk = res(xk);

% Compute objective and gradient at xk
fk = 0.5*(rk'*rk);  %avoid calling function func to not recompute res(xk)
gk = jk'*rk;        %avoid calling function grad_func to not recompute res(xk) and jac(xk)

num_of_iterations = 0;
pcg_iterations = [];

% % Exit imediatly if convergence criterion is satisfied
% if norm(gk,'inf')<tol
%     residual = rk;
%     exit_flag = 3;
%     fprintf('Critical point at ||g|| =  %.2g\n',norm(gk));
%     fprintf('LM time: %.2f\n', toc);
%     return;
% end

num_of_failed_steps = 0;

if verbose
report_diagonostic_info(0, norm(rk),norm(gk), lambda);
end

tic;

while true

    %rk,jk,fk and gk have already been computed
    
    % Matrix-free operators ( assemble A and b only once outside pcg)
    Afun  = @(x) jk' * (jk*x) + lambda * rt_r(x,jk);
    b = - gk;

    % Solve Afun * p = b (use pcg with optional preconditioner M)
    if isempty(M)
        [pk,flag,relres,iter] = pcg(Afun, b, tol, numel(x0));
        pcg_iterations(num_of_iterations+1) = iter;
    else
        P = @(v) M(v,lambda);
        [pk,flag,relres,iter] = pcg(Afun, b, tol, numel(x0), P);
        pcg_iterations(num_of_iterations+1) = iter;
    end

    %Line search
    % sc.p = xk;
    % sc.f = fk;
    % sc.g = gk;
    % [alpha,sp] = step_size_strong_wolfe(f,grad,sc,pk,1.0,params);

    %Try no line search (this is faster and gives better results for the
    %example considered)
    alpha = 1;
    
    % Check trial solution
    x_trial = xk + alpha*pk;

    r_trial = res(x_trial);
    f_trial = 0.5 * (r_trial'*r_trial);
    
    % Actual reduction in cost function
    actual = fk - f_trial;
    
    % Predicted reduction using quadratic model
    pred = 0.5*pk'*(lambda*pk-gk);
    
    % Safety: avoid dividing by zero / negative predicted reduction
    if abs(pred)<= 1e-15
        rho = -inf;   % treat as unsuccessful
    else
        rho = actual / pred;
    end

    % Levenberg–Marquardt update of lambda
    if rho > 0        % step is successful
        
        %Update lambda
        lambda = lambda * max(1/3, 1 - (2*rho - 1)^3);
        lambda = max(lambda, 1e-12);

        xk = x_trial;

        % Jacobian and residual at xk
        jk = jac(xk);
        rk = r_trial; %no need to recompute rk = res(xk);
        
        % Compute objective and gradient at xk
        fk = 0.5*(rk'*rk);  %avoid calling function func to not recompute res(xk)
        gk = jk'*rk;        %avoid calling function grad_func to not recompute res(xk) and jac(xk)

    else              % step is rejected
        lambda = lambda * nu;
        nu = 2*nu;
        num_of_failed_steps = num_of_failed_steps+1;
    end
    
    % Update counter 
    num_of_iterations = num_of_iterations + 1;
    
    if verbose
        report_diagonostic_info(num_of_iterations, norm(rk),norm(gk), lambda);
    end
    
    % Stopping Criterion (1) : maximum number of iterations
    if num_of_iterations >= max_iterations
        residual = rk;
        exit_flag = 1;
        fprintf('LM time: %.2f\n', toc);
        fprintf('LM exceeded max # of iterations: %i\n', max_iterations);
        return;
    end
    
    % Stopping Criterion (2) : maximum number of failed steps
    if num_of_failed_steps>= max_num_of_failed_steps
        residual = rk;
        exit_flag = 2;
        fprintf('LM time: %.2f\n', toc);
        fprintf('LM exceeded max # of FAILED steps: %i\n', max_num_of_failed_steps);
        return;
    end

    % Stopping criterion (3): Gradient of cost function is small
    if norm(gk,'inf')<tol
        residual = rk;
        exit_flag = 3;
        fprintf('Critical point at ||g|| =  %.2g\n',norm(gk));
        fprintf('LM time: %.2f\n', toc);
        return;
    end

    % Stopping criterion (4): Minimum number of steps, but no change in solution
    if num_of_iterations>min_num_of_steps && norm(pk)<=tol*(norm(xk)+tol)
        residual = rk;
        exit_flag = 4;
        fprintf('No change in solution ||pk|| = %.2g\n',norm(pk));
        fprintf('LM time: %.2f\n', toc);
        return;
    end

    % Stopping criterion (5): Cost function no longer decreasing
    if abs(fk-f_trial)<=tol && num_of_iterations>min_num_of_steps
        residual = rk;
        exit_flag = 5;
        fprintf('No change in cost function\n');
        fprintf('LM time: %.2f\n', toc);
        return
    end

end

end



% ------------------ diagnostics ---------------------
function report_diagonostic_info(iteration, residual_norm,grad_norm, lambda)
if nargin < 3
    fprintf('LM | it %i | r=%.5g| g=%.5g\n', iteration, residual_norm,grad_norm);
else
    fprintf('LM | it %i | r=%.5g | g=%.5g | lambda=%.2e\n', ...
        iteration, residual_norm,grad_norm, lambda);
end
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