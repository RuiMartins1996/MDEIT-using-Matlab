%Compute r matrix that measures distance between magnetometers and elements (depends only on geometry)
function [Rx,Ry,Rz] = computeRmatrices(fmdl,elementCentroids,sensorLocations)
    numElements = size(fmdl.elems,1);
    numSensors = size(sensorLocations,1);

    Rx = zeros(numSensors,numElements);
    Ry = zeros(numSensors,numElements);
    Rz = zeros(numSensors,numElements);
    
    for m = 1:numSensors
            sensorCenter = [sensorLocations(m,1),sensorLocations(m,2),sensorLocations(m,3)];
            for k = 1:numElements
                elementCentroid = elementCentroids(k,:);

                %Vector from centroid of element to center of sensor
                r = sensorCenter - elementCentroid;
                
                norm3r = norm(r)^3;

                Rx(m,k) = r(1)/norm3r;
                Ry(m,k) = r(2)/norm3r;
                Rz(m,k) = r(3)/norm3r;
            end
    end

end
