classdef testGNandLM < matlab.unittest.TestCase
    % Tests for Gauss–Newton and Levenberg–Marquardt solvers
    properties
        test_parameters;
    end

    methods(TestClassSetup)
        function setup(test_case)
            
            % Define test parameters
            
            A = [1 2; 3 5;-1 3];
            b = [3;1;0];

            res = @(x) A*x-b;
            
            jac = @(x) A;


            test_case.test_parameters.res = res;
            test_case.test_parameters.jac = jac;
            test_case.test_parameters.sol = inv(A'*A)*A'*b;
        end
    end

    methods (Test)
        
        function test_gauss_newton(test_case)

            res = test_case.test_parameters.res;
            jac = test_case.test_parameters.jac;

            % Nonlinear test problem
            
            x0 = [0;0];
            tol = 1e-8;
            maxit = 50;
            RtR = @(x,J) eye(numel(x));       % Identity Tikhonov
            lambda = 1e-5;                % GN (no damping)
            M = [];
            
            [xk, it, allx] = gn_solve(res, jac, ...
                x0, tol, maxit, RtR, lambda, M, false);
            
            sol = test_case.test_parameters.sol;
            fprintf('GN solution: [%.1f,%.1f]\n',xk(1),xk(2));
            fprintf('Linsolve solution: [%.1f,%.1f]\n',sol(1),sol(2));
            rk = res(xk);
            fprintf('Residual: [%.1f,%.1f,%.1f]\n',rk(1),rk(2),rk(3));
            
            test_case.verifyLessThan(norm(xk - test_case.test_parameters.sol(:)), numel(xk)*tol);
        end
        
        
        function test_levenberg_marquadt(test_case)

            res = test_case.test_parameters.res;
            jac = test_case.test_parameters.jac;

            % LM on same nonlinear problem
            
            x0 = [0;0];
            tol = 1e-8;
            maxit = 50;
            RtR = @(x,J) eye(numel(x));       % Identity Tikhonov
            lambda0 = 1e-2;
            M = [];

            [xk, it, allx] = lm_solve(res, jac, ...
                x0, tol, maxit, RtR, lambda0, M, false);

            sol = test_case.test_parameters.sol;
            fprintf('LM solution: [%.1f,%.1f]\n',xk(1),xk(2));
            fprintf('Linsolve solution: [%.1f,%.1f]\n',sol(1),sol(2));
            rk = res(xk);
            fprintf('Residual: [%.1f,%.1f,%.1f]\n',rk(1),rk(2),rk(3));


            test_case.verifyLessThan(norm(xk - test_case.test_parameters.sol(:)), numel(xk)*tol);

        end
        
    end
end
