function [xk,numOfIterations,allxk] = gaussNewton(res,jac,x0,tol,maxIterations,R,lambda)

% Pre - conditioner
if nargin <6
    R = speye(length(x0));
    lambda = 0.1;
end
   
xk = x0;

allxk = nan(length(x0),maxIterations+1);
allxk(:,1) = x0;

r0 = res(x0);

normR0 = norm(r0,2);

numOfIterations = 0;

k= 0;
while true
    numOfIterations = numOfIterations+1;
    
    % Compute resiudal
    r = res(xk);

    normRk = norm(r,2);

    fprintf('gn solver: iteration %i (r = %.5f)\n',numOfIterations,normRk)

    % Stopping criterion
    if normRk<tol
        break
    end
    
    % Stability check
    if normR0<normRk
        error('Bad step. Diverging!');
    end

    % Compute Jacobian
    J = jac(xk);
    
    % Construct NOSER prior
    R = sparse(1:length(x0),1:length(x0),sum(J.^2,1),length(x0),length(x0));
    
    %Compute search direction
    A = J'*J+lambda^2*(R'*R);
    b= J'*-r;
    pk = conjgrad(A,b,tol);
    
    %THERE IS A PROBLEM WITH THIS
    % SOLVER!!!! MATLAB DOING J'*J\(J'*-r) is way better, dunno why
    % pk = J'*J\(J'*-r);
    %But this will not work well with singular matrices
    %pk = (J'*J+lambda^2*(R'*R))\(J'*-r); 

    
    % Update solution
    xk = xk + pk;
    
    k=k+1;

    allxk(:,k+1) = xk;

    if k>maxIterations
        warning('Exceeded maximum number of iterations');
        break;
    end
end

end

