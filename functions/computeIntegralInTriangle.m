function integral = computeIntegralInTriangle(r,v1,v2,v3)

% Computes the Lkjm integral contribution to the magnetic field at point r

%% Check if inputs are 3D vectors
assert(isnumeric(v1) && isvector(v1) && numel(v1) == 3, ...
    'Input must be a 3D numeric vector.');
assert(isnumeric(v2) && isvector(v2) && numel(v2) == 3, ...
    'Input must be a 3D numeric vector.');
assert(isnumeric(v3) && isvector(v3) && numel(v3) == 3, ...
    'Input must be a 3D numeric vector.');
assert(isnumeric(r) && isvector(r) && numel(r) == 3, ...
    'Input must be a 3D numeric vector.');

% Check if points are colinear
tolerance = 1e-5;
assert(norm(cross(v2 - v1, v3-v1)) > tolerance,'Points are colinear!');

% For the first triangle, this gives the wrong result!
I(1) = computeIntegralInTriangleHelper(r,v1,v2,v3);

integral = I(1);

% x = [v1(1), v2(1), v3(1)];
% y = [v1(2), v2(2), v3(2)];
% z = [v1(3), v2(3), v3(3)];
% 
% % Plot filled triangle
% hold on
% fill3(x, y, z, 'magenta');  % 'magenta' is the fill color
% plot3(r(1),r(2),r(3),'b.','MarkerSize',10)
% axis equal;
% grid on;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% hold off



end


function integral = computeIntegralInTriangleHelper(r,v1,v2,v3)
%% Compute the auxiliary variables
D = [dot(v2-v1,v3-v1)/norm(v2-v1,2);...
    dot(v2-v1,v3-v2)/norm(v2-v1,2)];


% tri = [v1,v2,v3];
% fill3(tri(1,:),tri(2,:),tri(3,:),'g');

if abs(D(2))>1e-3 && abs(D(1))>1e-3

    A = [norm(v3-v1,2)^2;...
        norm(v3-v2,2)^2];

    B = [-2*dot(v3-v1,r-v1);...
        -2*dot(v3-v2,r-v2)];
    
    C = [norm(r-v1)^2;...
        norm(r-v2)^2];
    
    E = - [dot(v2-v1,r-v1)/norm(v2-v1);...
        dot(v2-v1,r-v2)/norm(v2-v1)];
    
    R = @(eta) sqrt(A.*eta^2+B.*eta+C);
    
    % These become huge numbers when the problematic elements are being
    % computed, why? It's because D(1) is very close to zero, but still
    % bigger than 1e-10. This happens because dot(v2-v1,v3-v1) =
    % -1.4618e-08, which means it's an almost perpendicular triangle. It
    % makes sense to relax to this condition abs(D(2))>1e-3 &&
    % abs(D(1))>1e-3, since even for a very i'll conditioned right
    % triangle. Assume there is an angle which is very close to 90 and is problematic.
    % Assuming that the triangles are not ill-conditioned, since
    % the meshing takes care of those, then the other angles will be
    % considerably different than 90. So it's okay to test the condition
    % for 1e-3 rather than for 1e-10. Maybe we should test directly
    % dot(v2-v1,v3-v1) and dot(v2-v1,v3-v2)???????
    
    a = A./(D.^2);

    b = (B.*D-2*A.*E)./(D.^2);

    c = (C.*D.^2-B.*D.*E+A.*E.^2)./(D.^2);

    
    % This expression has no numerical problems, while the above one does
    % d2 = c./(a-1) - (b.^2)./(4*(a-1).^2);
    d2 = (4.*c.*(a-1)-b.^2)./(4.*(a-1).^2);


    % To see if it is numerically zero, should check the numerator of this
    % expression:
    
    numerator_d2 = 4.*c.*(a-1)-b.^2;
    
    if abs(numerator_d2(1))<1e-14 d2(1) = 0; end
    if abs(numerator_d2(2))<1e-14 d2(2) = 0; end
    
    if numerator_d2(1)<0 d2(1) = 0; end
    if numerator_d2(2)<0 d2(2) = 0; end

    d = sqrt(d2);

    x = @(eta) D*eta+E;
    u = @(eta) x(eta)+b./(2*(a-1));

    %Compute the integral
    const = norm(cross(v2-v1,v3-v1),2)/norm(v2-v1,2);

    condition1 = R(1)+x(1)<1e-14 | R(0)+x(0)<1e-14;
    condition2 = abs(2*sqrt(a).*R(1)+2*a.*x(1)+b) < 1e-14 | abs(2*sqrt(a).*R(0)+2*a.*x(0)+b) < 1e-14;
    
    eta = 1;
    terms_at_1 = zeros(2,4);
    terms_at_1(:,1) = u(eta).*log(R(eta)+x(eta));
    terms_at_1(:,2) = -b./(2*sqrt(a).*(a-1)).*log(abs(2*sqrt(a).*R(eta)+2*a.*x(eta)+b));
    terms_at_1(:,3) = +d.*atan2(u(eta),d);
    terms_at_1(:,4) = -d.*atan2(2*d.*R(eta).*(a-1),(b.*x(eta)+2*c));
    
    terms_at_1(isnan(terms_at_1))=0;

    expr_at_1 = sum(terms_at_1,2);
    
    eta = 0;
    terms_at_0 = zeros(2,4);
    terms_at_0(:,1) = u(eta).*log(R(eta)+x(eta));
    terms_at_0(:,2) = -b./(2*sqrt(a).*(a-1)).*log(abs(2*sqrt(a).*R(eta)+2*a.*x(eta)+b));
    terms_at_0(:,3) = +d.*atan2(u(eta),d);
    terms_at_0(:,4) = -d.*atan2(2*d.*R(eta).*(a-1),(b.*x(eta)+2*c));

    terms_at_0(isnan(terms_at_0))=0;

    expr_at_0 = sum(terms_at_0,2);
    

    integral = const * sum([-1;1]./D*(1).*(expr_at_1-expr_at_0));
else
    integral = computeIntegralInTriangleHelper(r,v2,v3,v1);
end

end


function I = integrateOverTriangleCentroid(v1, v2, v3, rm)
% Ensure column vectors
v1 = v1(:); v2 = v2(:); v3 = v3(:); rm = rm(:);

% Compute the centroid of the triangle
centroid = (v1 + v2 + v3) / 3;

% Evaluate the integrand at the centroid
integrand = 1 / norm(centroid - rm);

% Compute the area of the triangle
normal = cross(v2 - v1, v3 - v1);
area = 0.5 * norm(normal);

% Approximate the integral as area * value at centroid
I = area * integrand;
end


function angles = triangle_angles_from_vertices(v1, v2, v3)
% TRIANGLE_ANGLES_FROM_VERTICES computes the interior angles of a triangle
% Input:
%   v1, v2, v3 - 1x2 or 1x3 vectors representing triangle vertices
% Output:
%   angles - 1x3 vector of interior angles (in radians), at [v1, v2, v3]

    % Side lengths opposite to each vertex
    a = norm(v2 - v3); % opposite to v1
    b = norm(v3 - v1); % opposite to v2
    c = norm(v1 - v2); % opposite to v3

    % Use Law of Cosines to find each angle
    angle1 = acos((b^2 + c^2 - a^2) / (2 * b * c)); % at v1
    angle2 = acos((c^2 + a^2 - b^2) / (2 * c * a)); % at v2
    angle3 = acos((a^2 + b^2 - c^2) / (2 * a * b)); % at v3

    angles = [angle1, angle2, angle3];
end