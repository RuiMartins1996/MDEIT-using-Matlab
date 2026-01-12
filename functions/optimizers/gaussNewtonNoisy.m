function [xk,numOfIterations,allxk] = ...
    gaussNewtonNoisy(res,jac,x0,tol,maxIterations,Wsqrt,RtR,lambda)

if not(isa(RtR, 'function_handle'))
    error('RtR should be a function handle')
end

if isempty(Wsqrt)
    Wsqrt = 1;
end


% We are going to adapt lambda according to Levenberg-Marquadt's method

% Initial damping parameter (according to page 25 of METHODS FOR NON-LINEAR
% LEAST SQUARES PROBLEMS Madsen, Nielsen, Tingleff)

if nargin < 8 || isempty(lambda) 
    lambda = 1;
end

    
lambda_inc = 10;    % factor to increase if step rejected
lambda_dec = 0.1;   % factor to decrease if step accepted

lambda_min = 1e-12; % safeguard
lambda_max = 1e+6;  % safeguard

% rho_min = 0.25;
% rho_max = 0.75;

% Constant regularization!
rho_min = -Inf;
rho_max = Inf;

xk = x0;

allxk = nan(length(x0),maxIterations+1);
allxk(:,1) = x0;

numOfIterations = 0;

k= 0;

% Compute residual
rk = res(xk);
% Prewhiten residual
rk = Wsqrt*rk;
% The objective function
phi_old = 0.5*norm(rk,2)^2;

fprintf('LM | it %i | r_prewhitened =%.5f \n',0,norm(rk,2));

% Discrepancy principle
tau = 1.1;

% delta = || y^\delta - y^\true || = noiseStdDev*sqrt(length(rk));
% Don't actually need delta when residual is prewhitened

numOfConsecutiveFailedSteps = 0;

while true
    
    % Compute Jacobian
    Jk = jac(xk);

    % Prewhiten Jacobian
    Jk = Wsqrt*Jk;
    
    % lhs matrix 
    A = Jk'*Jk+lambda*RtR(Jk);
    b = Jk'*-rk;
       
    % Simple Jacobi preconditioner
    L = sqrt(diag(diag(A)));
    U = L;

    pk = pcg(A,b,tol,numel(x0),L,U);

    % Trial step
    x_trial = xk + pk;
    rk_trial = res(x_trial);
    rk_trial_w = Wsqrt * rk_trial;
    phi_new = 0.5*norm(rk_trial_w,2)^2;
    
    % Predicted decrease from quadratic model
    pred = 0.5*pk'*(lambda*pk-Jk'*rk);
    
    % Compute rho
    rho = (phi_old - phi_new)/pred;
    
    if rho < rho_min
        % Reject step and increase lambda
        lambda = min(lambda*lambda_inc, lambda_max);

        % Stop if the consecutive number of failed steps is 10
        numOfConsecutiveFailedSteps =  numOfConsecutiveFailedSteps+1;
        if numOfConsecutiveFailedSteps >= 10
            error('Number of consecutive failed steps is 10. Exiting');
        end
    else 
        % Reset failed steps
        numOfConsecutiveFailedSteps = 0;
        
        % Update counter
        k=k+1;
        numOfIterations = numOfIterations+1;

        % Report diagonostic info
        fprintf('LM | it %i | r_prewhitened=%.5f | lambda=%.1e\n',...
            numOfIterations,norm(rk_trial_w,2),lambda);

        % % Stopping criterion (discrepancy principle applied to the
        % % pre-whitened residual)
        % if ...
        %     norm(rk_trial_w,2) <= tau*sqrt(length(rk))... % discrepancy principle
        %     % abs((norm(rk_trial_w,2)-norm(rk,2)))/norm(rk,2)<tol      % Residual change between iterations small
        %     fprintf('Converged at iteration %i with prewhitened residual %.5e\n', ...
        %         numOfIterations, norm(rk_trial_w,2));
        %     break;
        % elseif k>=maxIterations
        %     fprintf('LM exceeded max # of iterations: %i\n',maxIterations);
        %     break
        % end

        
        % Accept step
        xk = x_trial;

        % Save state in matrix
        allxk(:,k+1) = xk;

        if k>=maxIterations
            fprintf('LM exceeded max # of iterations: %i\n',maxIterations);
            break
        end

        %Update objective function
        phi_old = phi_new;

        % Update prewhitened residual
        rk =rk_trial_w;

        % Step very successful, reduce damping and accept
        if rho > rho_max
            % Decrease lambda (more Gauss-Newton-like)
            lambda = max(lambda*lambda_dec, lambda_min);
        else % rho within acceptable range, keep lambda constant and accept

        end
    end
end

end

