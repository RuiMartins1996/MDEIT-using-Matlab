function [xk,numOfIterations,allxk] = gaussNewtonLineSearch(res,jac,x0,tol,maxIterations)

xk = x0;

allxk = nan(length(x0),maxIterations+1);
allxk(:,1) = x0;

r0 = res(x0);
normR0 = norm(r0,2);

numOfIterations = 0;

% This should be updated later. The rank does not need to be computed like
% this. We know the rank as a function of the number of measurements
Jk = jac(x0);
rankJ = rank(Jk); % Rank of the Jacobian matrix

k= 0;
while true
    numOfIterations = numOfIterations+1;

    % Compute resiudal
    rk = res(xk);

    normRk = norm(rk,2);

    fprintf('gn | it %i | (||r|| = %.5f)\n',numOfIterations,normRk)

    % Stopping criterion
    if normRk<tol
        break
    end
    
    %Stability check
    if normR0<normRk
        warning('Bad step. Diverging!');
    end
    
    % Compute Jacobian
    Jk = jac(xk);
    
    % Since J^T*J is symmetric, then in the eigenvalue decomposition V*D*V^T,
    % it's true that V^T*V = 1

    % The singular value decomposition of J is U S V^T, where V is the same
    % as above. 

    [U,S,V] = svds(Jk,rankJ,'largest');
    
    %this comes from working out the GN normal equations with the SVD for J 
    pk = -V*((U'*rk)./diag(S)); 
    
    % Line search over this direction 

    % Objective and gradient
    func = @(x) 0.5 * (res(x)' * res(x));
    grad = @(x) jac(x)' * res(x);

    fk =  0.5 *(rk'*rk);
    gk = Jk'*rk;

    alphak = strongWolfeConditions(func,grad,xk,fk,gk,pk);
    
    % Update solution
    xk = xk + alphak*pk;
    
    k=k+1;

    allxk(:,k+1) = xk;

    if k>maxIterations
        normRk = norm(res(xk),2);
        fprintf('gn | iteration %i | (r = %.5f)\n',numOfIterations,normRk)
        fprintf('gn exceeded max # of iterations: %i\n',maxIterations);
        break;
    end
end

end



%% FUNCTION: Strong Wolfe Conditions
function alphastar = strongWolfeConditions(func,grad,x0,f0,g0,p)
% function [alpha] = strong_wolfe(func,x0,f0,g0,p)

% Parameters
c1 = 1e-4;
c2 = 0.9;
alphaMax = 2.5;
maxIterations = 20;

fim1 = f0;
alphai = 1.0;
alphaim1 = 0.0;
dphi0 = g0'*p;

% search for alpha that satisfies strong-Wolfe conditions
i=1;
while true

    xi = x0 + alphai*p;

    %Evaluate cost funtion at xi
    fi = func(xi);

    if (fi > f0 + c1*alphai*dphi0) || ( (fi >= fim1) && i>1 )
        gi = grad(xi);
        alphastar = zoom(func,grad,x0,f0,g0,p,alphaim1,alphai);
        break;
    end

    gi = grad(xi);
    dphi = gi'*p;
    if abs(dphi) <= -c2*dphi0
        alphastar = alphai;
        break
    end

    if dphi >=0
        alphastar = zoom(func,grad,x0,f0,g0,p,alphai,alphaim1);
        break;
    end

    alphai = alphai + 0.8*(alphaMax-alphai);

    if (i > maxIterations)
        alphastar = alphai;
        warning('Line search exceeded maximum number of iterations!');
        break;
    end

    i = i+1;
end

end

%% FUNCTION: Zoom
function alphastar = zoom(func,grad,xi,fi,gi,p,alphaLo,alphaHi)

% Parameters
c1 = 1e-4;
c2 = 0.9;
maxIterations = 20;

dphi0 = gi'*p;

j = 1;
while true
    alphaj = 0.5*(alphaLo + alphaHi);
    
    xj = xi + alphaj*p;
    xlo = xi + alphaLo*p;

    fj = func(xj);
    flo = func(xlo);
    
    if fj > fi + c1*alphaj*dphi0 || fj >= flo
        alphaHi = alphaj;
    else
        gj = grad(xj);
        dphj = gj'*p;
        if abs(dphj) <= -c2*dphi0 
            alphastar = alphaj;
            break;
        end
        if dphj * (alphaHi-alphaLo) >= 0
            alphaHi = alphaLo;
        end
        alphaLo = alphaj;
    end
    j = j+1;

    if (j > maxIterations)
        alphastar = alphaj;
        break;
    end
end

end