function elec = assemble_electrodes(ctx, theta_minus, width, z_contact)
%ASSEMBLE_ELECTRODES  Build the M electrode arc integrals + gap geometry
%for a design vector theta_minus (PLAN_implementation.md sec 3.2-3.3).
%
%   elec = assemble_electrodes(ctx, theta_minus, width, z_contact)
%
% theta_minus: [M x 1], the left endpoints in the FIXED CYCLIC ORDER
%   1,2,...,M around the circle (electrode identity = index, not sorted
%   theta value -- the repulsion term in the cost keeps them from
%   crossing during optimization). Not reduced mod 2pi -- keep it as a
%   monotonically increasing running angle across the M electrodes so
%   ordering/gaps are unambiguous; only reduced mod 2pi inside
%   electrode_arc_integrals when locating mesh edges.
% width: fixed arc-length electrode width (scalar, same for all M).
% z_contact: scalar or [M x 1] contact impedances.
%
% Returns a 1xM struct array elec(m) with fields:
%   theta_minus, theta_plus, dtp_dtm   (scalars)
%   z                                   contact impedance
%   Bm, cm, Elen                        electrode arc integrals
%   phi_minus, phi_plus, speed_minus, speed_plus   (for FTC derivatives)
%   gap_next    arc length from this electrode's right end to the NEXT
%               electrode's left end (cyclic, electrode M wraps to
%               electrode 1 shifted by +2pi)

M = numel(theta_minus);
shape = ctx.shape;
if isscalar(z_contact), z_contact = z_contact*ones(M,1); end

[theta_plus, dtp_dtm] = solve_right_endpoint(theta_minus, width, shape);

elec = struct('theta_minus',cell(1,M),'theta_plus',cell(1,M), ...
    'dtp_dtm',cell(1,M),'z',cell(1,M),'Bm',cell(1,M),'cm',cell(1,M), ...
    'Elen',cell(1,M),'phi_minus',cell(1,M),'phi_plus',cell(1,M), ...
    'speed_minus',cell(1,M),'speed_plus',cell(1,M),'gap_next',cell(1,M));

for m = 1:M
    E = electrode_arc_integrals(ctx, theta_minus(m), theta_plus(m));
    elec(m).theta_minus = theta_minus(m);
    elec(m).theta_plus  = theta_plus(m);
    elec(m).dtp_dtm     = dtp_dtm(m);
    elec(m).z    = z_contact(m);
    elec(m).Bm   = E.Bm;
    elec(m).cm   = E.cm;
    elec(m).Elen = E.Elen;
    elec(m).phi_minus = E.phi_minus;
    elec(m).phi_plus  = E.phi_plus;
    elec(m).speed_minus = E.speed_minus;
    elec(m).speed_plus  = E.speed_plus;
end

g = electrode_gaps(theta_minus, width, shape);
for m = 1:M
    elec(m).gap_next = g(m);
end
end
