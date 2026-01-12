function [xk,numOfIterations] = steepestDescent(f,grad,x0,tol,maxIterations)

xk = x0;

numOfIterations = 0;

colors = ['r','g','b','y','k','c','m'];

color = colors(randi(length(colors)));
for k = 1:maxIterations

    numOfIterations = numOfIterations+1;
    plot(xk(1),xk(2),'o','MarkerSize',5,'Color',color);

    %Compute gradient at xk
    dfx = grad(xk);
    if size(dfx,2)>size(dfx,1)
        dfx = dfx';
    end

    %Stopping criterion
    if f(xk)<tol
        break
    end

    %Search direction
    pk = -dfx;

    quiver(xk(1),xk(2),pk(1),pk(2),'k-')
    pause(1e-10)

    alphak = backtrackingLineSearch(f,grad,pk,xk);

    xk = xk + alphak*pk;
end

if k>=maxIterations
    warning('Steepest Descent with FD exceeded maximum number of iterations');
end

end

