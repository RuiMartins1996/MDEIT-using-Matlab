function [L,n] = computeIntegralAndNormalsInTetrahedron(vertices,rm)

assert(all(size(vertices)==[4,3]),'Input does not have the correct size');
assert(size(rm,1) == 1, 'Input rm must be a line vector');

faces = [1   2   3;...
        1   2   4;...
        1   3   4;...
        2   3   4];

n = computeOutwardUnitNormals(vertices,faces);

% Compute the integral on each triangle
L = zeros(1,4);

for i = 1:size(faces,1)

     v1 = vertices(faces(i,1),:);
     v2 = vertices(faces(i,2),:);
     v3 = vertices(faces(i,3),:);

     L(i) = computeIntegralInTriangle(rm,v1,v2,v3);
end

end


function n = computeOutwardUnitNormals(vertices,faces)

% Compute tetrahedron centroid (volume center)
tetra_centroid = mean(vertices, 1);

n = zeros(3,4);

% Compute and plot unit normal vectors
for i = 1:size(faces,1)
    idx = faces(i,:);
    p1 = vertices(idx(1),:);
    p2 = vertices(idx(2),:);
    p3 = vertices(idx(3),:);
    
    % Two edge vectors
    v1 = p2 - p1;
    v2 = p3 - p1;

    % Face normal using cross product
    ni = cross(v1, v2);
    ni = ni / norm(ni); % Normalize

    % Face centroid
    face_centroid = (p1 + p2 + p3) / 3;

    % Vector from face centroid to tetrahedron centroid
    to_center = tetra_centroid - face_centroid;

    % Flip normal if it's pointing inward
    if dot(ni, to_center) > 0
        ni = -ni;
    end

    n(:,i) = ni;
    
end

return

end