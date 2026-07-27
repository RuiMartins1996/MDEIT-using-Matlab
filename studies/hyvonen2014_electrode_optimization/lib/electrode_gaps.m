function g = electrode_gaps(theta_minus, width, shape)
%ELECTRODE_GAPS  Cyclic arc-length gaps g_m between electrode m's right
%end and electrode m+1's left end (electrode M wraps to electrode 1).
%Pure geometry -- no mesh/FEM solve -- so it is cheap enough to call
%repeatedly inside the Algorithm-1 line-search bisection
%(steepest_descent_algorithm1.m) and inside assemble_electrodes.m.
%
%   g = electrode_gaps(theta_minus, width, shape)   % [M x 1]

M = numel(theta_minus);
theta_plus = solve_right_endpoint(theta_minus, width, shape);

g = zeros(M,1);
for m = 1:M
    mnext = mod(m,M)+1;
    if mnext == 1
        target = theta_minus(1) + 2*pi;
    else
        target = theta_minus(mnext);
    end
    g(m) = arc_length_between(theta_plus(m), target, shape);
end
end
