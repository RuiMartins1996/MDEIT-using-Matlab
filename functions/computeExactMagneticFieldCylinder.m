function B = computeExactMagneticFieldCylinder(rQ,R,L,J0)

global mu0;
% Assume that the cylinder is centered at (0,0,0) with height spanning
% [-L/2,L/2]

% rQ    : query point
% R     : cylinder radius
% L     : cylinder height
% J0    : DC current density going trough the cylinder

zC = rQ(3);
RC = sqrt(rQ(1)^2+rQ(2)^2);

if zC<0
    rQNew = rQ;
    rQNew(3) = -rQ(3);
    B = computeExactMagneticFieldCylinder(rQNew,R,L,J0);
    return;
end

AbsTol = 1e-14;
RelTol = 1e-14;


%% Define the sigma+- functions

hPlus = sqrt((L/2+zC)^2);
hMinus = sqrt((L/2-zC)^2);

integrandPlus = @(r,theta) ...
    (hPlus ./ (hPlus^2 + r.^2 + RC^2 + 2*r*RC.*cos(theta))) .* ...
    ((hPlus^2 + r.^2 + RC*r.*cos(theta)) ./ (hPlus^2 + r.^2.*(sin(theta)).^2));

integrandMinus = @(r,theta) ...
    (hMinus ./ (hMinus^2 + r.^2 + RC^2 + 2*r*RC.*cos(theta))) .* ...
    ((hMinus^2 + r.^2 + RC*r.*cos(theta)) ./ (hMinus^2 + r.^2.*(sin(theta)).^2));

sigmaPlus = @(r) 1/pi*integral(@(theta) integrandPlus(r,theta), 0, pi,'ArrayValued',true,'RelTol',RelTol,'AbsTol',AbsTol);
sigmaMinus = @(r) 1/pi * integral(@(theta) integrandMinus(r,theta),0,pi,'ArrayValued',true,'RelTol',RelTol,'AbsTol',AbsTol);

%% Compute the magnetic field intensity

BC = nan;

% First region of interest
if zC<=L/2 && zC>=0

    % Case 1: RC>R (equation 19)
    if RC>R
        BC = mu0*J0/(2*RC)*integral(@(r) r.*(sigmaPlus(r)+sigmaMinus(r)), 0, R,'RelTol',RelTol,'AbsTol',AbsTol);
    end

    % Case 2: RC<R (equation 23)
    if RC<R
        BC = mu0*J0/(2*RC)*(...
            -2*integral(@(r) r, RC, R)+...
            integral(@(r) r.*(sigmaPlus(r)+sigmaMinus(r)), 0, R,'RelTol',RelTol,'AbsTol',AbsTol)...
            );
    end
end

% Second region of interest
if zC>=L/2
    BC = mu0*J0/(2*RC)*(...
        integral(@(r) r.*(sigmaPlus(r)-sigmaMinus(r)), 0, R,'RelTol',RelTol,'AbsTol',AbsTol)...
        );
end

%% Output the magnetic field vector

% B is tangential to the C circle

t = [-rQ(2) rQ(1) 0];
t = t/norm(t);

B = BC*t;

end