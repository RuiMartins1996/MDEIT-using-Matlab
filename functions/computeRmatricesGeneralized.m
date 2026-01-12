function [Rx,Ry,Rz] = computeRmatricesGeneralized(fmdl,sensorLocations)

numElements = size(fmdl.elems,1);
numSensors = size(sensorLocations,1);

Rx = zeros(numSensors,numElements);
Ry = zeros(numSensors,numElements);
Rz = zeros(numSensors,numElements);


for m = 1:numSensors
    sensorCenter = [sensorLocations(m,1),sensorLocations(m,2),sensorLocations(m,3)];
    
    for k = 1:numElements
        
        element_ids = fmdl.elems(k,:);

        vertices = fmdl.nodes(element_ids,:);

        [L,n] = computeIntegralAndNormalsInTetrahedron(vertices,sensorCenter);
        
        S = zeros(3,1);
        for j = 1:4
            S = S + L(j)*n(:,j);
        end

        Rx(m,k) = dot(S,[1,0,0]);
        Ry(m,k) = dot(S,[0,1,0]);
        Rz(m,k) = dot(S,[0,0,1]);
    end

end

end