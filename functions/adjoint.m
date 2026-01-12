%FUNCTION: This will solve 
function df = adjoint(x,p,A,lambdatimesdAdp,lambdatimesdbdp,dfdx,dfdp,tol)


if isa(A, 'function_handle')
    b = -dfdx(x,p);
    if size(b,2)>size(b,1)
        b=b';
    end

    % Simple Jacobi preconditioner
    L = sqrt(diag(diag(A(p))));
    U = L;
    
    [lambda,flag,relres] = pcg(A(p),b,tol,numel(p),L,U);
else

    % Simple Jacobi preconditioner
    L = sqrt(diag(diag(A)));
    U = L;

    [lambda,flag,relres] = pcg(A,-dfdx(x,p),tol,numel(p),L,U);
end

if size(lambda,1)>size(lambda,2)
    lambda = lambda';
end

%df1 = tensorprod(lambdatimesdAdp(lambda,p),x,2,1);

df1 = x'*lambdatimesdAdp(lambda,p);

df2 = -lambdatimesdbdp(lambda,p);

b = dfdp(x,p);

if size(b,1)>size(b,2)
    b = b';
end

df = df1+df2+b;
% df = lambda*(tensorprod(dAdp(p),x,2,1)-dbdp(p))+dfdp(x,p)';

end

