function [theta_plus, dthetaplus_dthetaminus] = solve_right_endpoint(theta_minus, width, shape)
%SOLVE_RIGHT_ENDPOINT  Electrode right endpoint for a fixed arc-length width.
%
%   [theta_plus, dtp_dtm] = solve_right_endpoint(theta_minus, width, shape)
%
% Solves  int_{theta_minus}^{theta_plus} |gammadot(t)| dt = width  for
% theta_plus (theta_plus > theta_minus, arbitrary real, NOT wrapped to
% [0,2pi) -- callers keep the running/unwrapped value so that gaps and
% ordering between consecutive electrodes are unambiguous).
%
% dtp_dtm = dtheta_plus/dtheta_minus = |gammadot(theta_minus)| / |gammadot(theta_plus)|
% via (4.2) in the paper / eq (4.2) in PLAN_implementation.md.
%
% Fast path for the disk (|gammadot| = r, constant): theta_plus is
% analytic and the derivative is exactly 1. The general branch (bisection
% on the arc-length integral) is provided for future non-circular shapes
% (Case 3) but is not exercised/validated by Cases 1-2.

theta_minus = theta_minus(:);
n = numel(theta_minus);

switch shape.type
    case 'disk'
        r = shape.radius;
        theta_plus = theta_minus + width./r;
        dthetaplus_dthetaminus = ones(n,1);
    otherwise
        theta_plus = zeros(n,1);
        dthetaplus_dthetaminus = zeros(n,1);
        for k = 1:n
            f = @(tp) arc_length_between(theta_minus(k), tp, shape) - width;
            % bracket and bisect
            lo = theta_minus(k); hi = theta_minus(k) + 4*pi;
            while f(hi) < 0
                hi = hi + 2*pi;
            end
            for iter = 1:80
                mid = 0.5*(lo+hi);
                if f(mid) > 0, hi = mid; else, lo = mid; end
            end
            tp = 0.5*(lo+hi);
            theta_plus(k) = tp;
            [~,~,s_minus] = boundary_curve(theta_minus(k), shape);
            [~,~,s_plus]  = boundary_curve(tp, shape);
            dthetaplus_dthetaminus(k) = s_minus/s_plus;
        end
end
end
