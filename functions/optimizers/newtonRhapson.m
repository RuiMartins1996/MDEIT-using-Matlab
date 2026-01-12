function [xk,numOfIterations] = newtonRhapson(res,jac,x0,tol,maxIterations)
   
xk = x0;

numOfIterations = 0;

for k = 1:maxIterations

    plot(xk(1),xk(2),'o','MarkerSize',10);

    numOfIterations = numOfIterations+1;
    
    % Compute Jacobian
    J = jac(xk);

    % Compute resiudal
    r = res(xk);

    % Stopping criterion
    if norm(r,inf)<tol
        break
    end

    %Compute search direction
    pk = J\-r;

    quiver(xk(1),xk(2),pk(1),pk(2),'k-')


    % Update solution
    xk = xk + pk;
end

end

