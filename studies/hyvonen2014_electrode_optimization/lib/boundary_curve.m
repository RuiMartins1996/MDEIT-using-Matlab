function [g, gdot, speed] = boundary_curve(theta, shape)
%BOUNDARY_CURVE  Star-shaped boundary parametrization gamma(theta).
%
%   [g, gdot, speed] = boundary_curve(theta, shape)
%
% theta may be a column vector of parameter values. shape is a struct with
% shape.type in {'disk'} (only 'disk' is implemented; this is the only
% shape used by Cases 1-2 of the plan). g is [numel(theta) x 2], gdot is
% the derivative dgamma/dtheta (same size), speed = |gdot| (column).
%
% Kept as a separate, generic function (rather than hard-coding cos/sin
% everywhere) so a future session can add 'ellipse'/'peanut' shapes for
% Case 3 without touching any other file -- every other routine in lib/
% only calls boundary_curve and solve_right_endpoint, never cos/sin
% directly.

theta = theta(:);

switch shape.type
    case 'disk'
        r = shape.radius;
        g    = r*[cos(theta), sin(theta)];
        gdot = r*[-sin(theta), cos(theta)];
        speed = r*ones(numel(theta),1);
    otherwise
        error('boundary_curve:unsupported','shape.type "%s" not implemented',shape.type);
end
end
