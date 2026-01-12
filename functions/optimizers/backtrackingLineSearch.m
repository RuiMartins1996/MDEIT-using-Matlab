function alpha = backtrackingLineSearch(f,grad,p,x)
tau = 0.5;
c1 = 1e-4;
alpha = 1;

dphi0 = dot(grad(x),p);
fx = f(x);

maxIterations = 20;

for j = 1:maxIterations
    
    % t = -c1*dphi0;
    if f(x+alpha*p) > fx + alpha*c1*dphi0
        alpha = tau*alpha;
    else
        break; 
    end
end

if j == maxIterations
    warning('Maximum number of iterations exceeded in backtrackingLineSearch');
end

return;
end