function [xk1,numOfIterations] = BFGS(f,grad,x0,tol,maxIterations)

N = length(x0);

%%  Cross-plataform way to determine free memory (AVOID USING VIRTUAL MEMORY)
runtime = java.lang.Runtime.getRuntime;

% Get the total memory available to Java (and MATLAB)

% totalMemory = runtime.totalMemory();  % in bytes
freeMemory = runtime.freeMemory();    % in bytes
freeMemoryGB = freeMemory / 1024^3;

% Calculate the size in gigabytes of inverse Hessian
sizeInGB = (N*N * 8) / (1024^3);

if sizeInGB > freeMemoryGB
    warning('on')
    warning('Size of inverse Hessian exceeds available RAM. Using virtual memory!') 
end

%% Initialize

%Initial guess for the inverse of the Hessian matrix
Hk = eye(N,N);

xk = x0;
fk = f(xk);
gk = grad(xk);

%Sanity check
checkInputs(xk,fk,gk);

%Display info
display = true;

if display
    fprintf(' iter        f(x)          |g(x)| \n')
    fprintf('-------------------------------------\n')
end

%% Iterate
numOfIterations = 1;

while numOfIterations<=maxIterations
    
    [bool,normgk] = checkOptimality(fk,gk,tol);
    if bool
        fprintf('%3d %16.8f %16.8f\n',numOfIterations,fk,normgk);
        break
    else
        fprintf('%3d %16.8f %16.8f\n',numOfIterations,fk,normgk);
    end

    %Compute search direction without solving linear system ( with inverse
    %of Bk)
    pk = -Hk*gk;

    % quiver(xk(1),xk(2),pk(1),pk(2),'k-')
    
    % Backtracking line search
    alphak = backtrackingLineSearch(f,grad,pk,xk);

    % alphak = strong_wolfe(f,grad,xk,fk,gk,pk);
    
    % Update inverse Hessian (Sherman-Morrison formula)
    sk = alphak*pk;

    xk1 = xk + sk;
    
    gk1 = grad(xk1);

    yk = gk1-gk;
    
    Hk = Hk + (sk'*yk+yk'*Hk*yk)*(sk*sk')/((sk'*yk)^2)-(Hk*yk*sk'+sk*yk'*Hk)/(sk'*yk);
    
    % % Check if solution is not advancing
    % opt = norm(sk,inf);
    % if opt<=tol
    %     break
    % end
     
    % Update solution
    xk = xk1;

    num_of_sensors = N/3;
    plot3(x0(1:num_of_sensors),x0(num_of_sensors+1:2*num_of_sensors),x0(2*num_of_sensors+1:3*num_of_sensors),'r.','MarkerSize',3)
    pause(1e-10)

    %Compute cost function at new xk for next iteration
    fk = f(xk);
    %Compute gradient at new xk for next iteration
    gk = grad(xk);

    numOfIterations = numOfIterations+1;
end

if numOfIterations > maxIterations
    warning('on')
    warning('Exceeded maximum number of iterations');
end

end

%% FUNCTION: CHECKINPUTS check if state vector and gradient are both column vectors
function checkInputs(x,f,g)
if not(isnumeric(f) && isscalar(f))
    error('Cost function f must return a number')
end

if size(x,2)>size(x,1)
    error('state vector x must return a column vector');
end

if size(g,2)>size(g,1)
    error('function grad must return a column vector');
end
end

%% FUNCTION: CHECKOPTIMALITY check if cost function is under specified tolerance

function [bool,normgk] = checkOptimality(f,g,tol)
    normgk = norm(g,inf);
    bool = f < tol || normgk<tol;
end

%% FUNCTIONS stolen from implementation of L-BFGS-B and changed

function [alpha] = strong_wolfe(func,grad,x0,f0,g0,p)
% function [alpha] = strong_wolfe(func,grad,x0,f0,g0,p)
% Compute a line search to satisfy the strong Wolfe conditions.
% Algorithm 3.5. Page 60. "Numerical Optimization". Nocedal & Wright.
% INPUTS:
%  func: objective function handle.
%  x0: [n,1] initial design vector.
%  f0: initial function evaluation.
%  g0: [n,1] initial objective gradient vector.
%  p: [n,1] search direction vector.
% OUTPUTS:
% alpha: search length

% initialize variables
c1 = 1e-4;
c2 = 0.9;
alpha_max = 2.5;
alpha_im1 = 0;
alpha_i = 1;
f_im1 = f0;
dphi0 = transpose(g0)*p;
i = 0;
max_iters = 20;

% search for alpha that satisfies strong-Wolfe conditions
while true
  
  x = x0 + alpha_i*p;
  f_i = feval(func, x);
  g_i = feval(grad, x);
  
  if (f_i > f0 + c1*dphi0) || ( (i > 1) && (f_i >= f_im1) )
    alpha = alpha_zoom(func,grad,x0,f0,g0,p,alpha_im1,alpha_i);
    break;
  end
  dphi = transpose(g_i)*p;
  if ( abs(dphi) <= -c2*dphi0 )
    alpha = alpha_i;
    break;
  end
  if ( dphi >= 0 )
    alpha = alpha_zoom(func,grad,x0,f0,g0,p,alpha_i,alpha_im1);
    break;
  end
  
  % update
  alpha_im1 = alpha_i;
  f_im1 = f_i;
  alpha_i = alpha_i + 0.8*(alpha_max-alpha_i);
  
  if (i > max_iters)
    alpha = alpha_i;
    break;
  end
  
  i = i+1;
  
end

end


function [alpha] = alpha_zoom(func,grad,x0,f0,g0,p,alpha_lo,alpha_hi)
% function [alpha] = alpha_zoom(func,x0,f0,g0,p,alpha_lo,alpha_hi)
% Algorithm 3.6, Page 61. "Numerical Optimization". Nocedal & Wright.
% INPUTS:
%  func: objective function handle.
%  x0: [n,1] initial design vector.
%  f0: initial objective value.
%  g0: [n,1] initial objective gradient vector.
%  p: [n,1] search direction vector.
%  alpha_lo: low water mark for alpha.
%  alpha_hi: high water mark for alpha.
% OUTPUTS:
%  alpha: zoomed in alpha.

% initialize variables
c1 = 1e-4;
c2 = 0.9;
i = 0;
max_iters = 20;
dphi0 = transpose(g0)*p;

while true
  alpha_i = 0.5*(alpha_lo + alpha_hi);
  alpha = alpha_i;
  x = x0 + alpha_i*p;
  
  f_i = feval(func, x);
  g_i = feval(grad, x);

  x_lo = x0 + alpha_lo*p;

  f_lo = feval(func, x_lo);
  
  if ( (f_i > f0 + c1*alpha_i*dphi0) || ( f_i >= f_lo) )
    alpha_hi = alpha_i;
  else
    dphi = transpose(g_i)*p;
    if ( ( abs(dphi) <= -c2*dphi0 ) )
      alpha = alpha_i;
      break;
    end
    if ( dphi * (alpha_hi-alpha_lo) >= 0 )
      alpha_hi = alpha_lo;
    end
    alpha_lo = alpha_i;
  end
  i = i+1;
  if (i > max_iters)
    alpha = alpha_i;
    break;
  end
end

end