function s = arc_length_between(t0, t1, shape)
%ARC_LENGTH_BETWEEN  int_{t0}^{t1} |gammadot(t)| dt  (t1 may be < t0; sign
%follows the integration direction). Shared by solve_right_endpoint.m and
%the electrode gap computation in assemble_electrodes.m.

switch shape.type
    case 'disk'
        s = shape.radius*(t1-t0);
    otherwise
        n_quad = 40;
        [xg,wg] = gauss_legendre_nodes(n_quad);
        t = t0 + 0.5*(t1-t0)*(xg+1);
        [~,~,speed] = boundary_curve(t, shape);
        s = 0.5*(t1-t0)*sum(wg.*speed);
end
end
