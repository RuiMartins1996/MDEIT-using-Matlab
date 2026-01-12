function img = computeGeometryMatrices(img,mu0)
fmdl = img.fwd_model;
sensorLocations = img.sensorLocations;
numOfSensors = size(sensorLocations,1);

% Compute and assign G and R
[Gx,Gy,Gz,V,elementCentroids] = computeGradientMatrix(fmdl);
% [Rx,Ry,Rz] = computeRmatrices(fmdl,elementCentroids,sensorLocations);

[Rx,Ry,Rz] = computeRmatricesGeneralized(fmdl,sensorLocations);
for k = 1:length(V)
   Rx(:,k) = Rx(:,k)./V(k);
   Ry(:,k) = Ry(:,k)./V(k);
   Rz(:,k) = Rz(:,k)./V(k);
end

img.G = struct('Gx',Gx,'Gy',Gy,'Gz',Gz);
img.R = struct('Rx',Rx,'Ry',Ry,'Rz',Rz);
img.elem_volume = V;
img.elem_centroids = elementCentroids;

% Assign mu0 (must go above: needed by computeGammaMatrices)
img.mu0 = mu0;

% Compute and assign dGamma
for m = 1:numOfSensors
    [dGammaX,dGammaY,dGammaZ] = ...
        computeGammaMatricesDerivatives(img,m);
    dGammaXcell{m} = dGammaX;
    dGammaYcell{m} = dGammaY;
    dGammaZcell{m} = dGammaZ;
end

img.dGammaXcell = dGammaXcell;
img.dGammaYcell = dGammaYcell;
img.dGammaZcell = dGammaZcell;
end