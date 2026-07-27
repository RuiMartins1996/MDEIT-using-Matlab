function check_gradient_fd(costgrad, theta0, n_check, h)
%CHECK_GRADIENT_FD  Central-difference gradient check (copied/adapted from
%studies/optimal_sensors_bayesian_approach/example_anomaly_circle_2d.m,
%the validation pattern this whole repo relies on for gradient-bug
%catching -- V4 in PLAN_implementation.md, the mandatory gate before any
%optimization run).
%
%   check_gradient_fd(costgrad, theta0, n_check, h)

n = numel(theta0);
[~,g] = costgrad(theta0);
idx = randperm(n,min(n_check,n));
for i = idx
    e = zeros(n,1); e(i) = 1;
    fp = costgrad(theta0 + h*e);
    fm = costgrad(theta0 - h*e);
    g_fd = (fp - fm)/(2*h);
    rel_err = abs(g(i)-g_fd)/max(abs(g_fd),eps);
    fprintf('  dpsi/dtheta(%i): analytic = %+.6e | FD = %+.6e | rel err = %.2e\n',...
        i,g(i),g_fd,rel_err);
end
end
