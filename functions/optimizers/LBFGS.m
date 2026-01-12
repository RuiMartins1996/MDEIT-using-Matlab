function [xkp1,k] = LBFGS(func,grad,x0,options)

% Parse options
[m,tol,maxIterations,display] = parseOptions(options);

n = length(x0);
xk = x0;

% Displacement vectors sk and gradient differences yk 
S = nan(n,m);
Y = nan(n,m);

k = 1;
bool = true;

if (display)
  fprintf(' iter        f(x)          optimality\n')
  fprintf('-------------------------------------\n')
end

while bool

    % plot(xk(1),xk(2),'k*')
    
    fk = func(xk);
    gk = grad(xk);
    
    [optimality,normgk] = checkOptimality(fk,gk,tol);

    if display
        fprintf('%3d %16.8f %16.8f\n',k,fk,normgk);
    end

    if optimality == true
        xkp1 = xk;
        break;
    end

    if k ==1
        gammak = 1;
    elseif k<=m
        gammak = S(:,k-1)'*Y(:,k-1)/(Y(:,k-1)'*Y(:,k-1));
    else
        gammak = S(:,end)'*Y(:,end)/(Y(:,end)'*Y(:,end));
    end

    H0k = gammak*speye(n,n);

    rk = twoLoopRecursion(gk,H0k,S,Y,min(m,k-1));
    pk = -rk;
    
    alphak = strongWolfeConditions(func,grad,xk,fk,gk,pk);
    
    xkp1 = xk + alphak*pk;
    sk = xkp1-xk;
    yk = grad(xkp1)-gk;
    
    if k>m
        S = [S(:,2:end) sk];
        Y = [Y(:,2:end) yk];
    else
        S(:,k) = sk;
        Y(:,k) = yk;
    end

    xk = xkp1;

    if k>=maxIterations
        warning('LBFGS exceeded maximum number of iterations');
        break;
    end

    k=k+1;
end


end

%% FUNCTION: Check optimality
function [bool,normgk] = checkOptimality(f,g,tol)
    normgk = norm(g,inf);
    bool = normgk<tol;
end

%% FUNCTION: Two loop recursion
function r = twoLoopRecursion(gk,H0k,S,Y,M)

alpha = zeros(1,M);
rho = zeros(1,M);

q = gk;

%k-1:-1:k-M

for i = M:-1:1
    rho(i) = 1/(S(:,i)'*Y(:,i));
    alpha(i) = rho(i)*S(:,i)'*q;
    q = q-alpha(i)*Y(:,i);
end

r = H0k*q;

for i = 1:M
    beta = rho(i)*Y(:,i)'*r;
    r = r + S(:,i)*(alpha(i)-beta);
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

%% FUNCTION: Parse options
function [m,tol,maxIterations,display] = parseOptions(options)

m = 10;
tol = 1.0e-5;
maxIterations = 20;
display = true;

if ( isfield(options, 'm') )
    m = options.m;
end
if ( isfield(options, 'tol') )
    tol = options.tol;
end
if ( isfield(options, 'maxIterations') )
    maxIterations = options.maxIterations;
end
if ( isfield(options, 'display') )
    warning('on');
    display = options.display;
end

end

