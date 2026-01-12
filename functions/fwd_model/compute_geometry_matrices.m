function fmdl = compute_geometry_matrices(fmdl,mu0)

% Compute and assign G and R ( calls the correct function for the 2d and 3d
% case)
[Gx,Gy,Gz,V,elementCentroids] = compute_gradient_matrix(fmdl);

[Rx,Ry,Rz] = compute_r_matrices(fmdl);

fmdl.G = struct('Gx',Gx,'Gy',Gy,'Gz',Gz);
fmdl.R = struct('Rx',Rx,'Ry',Ry,'Rz',Rz);
fmdl.elem_volume = V;
fmdl.elem_centroids = elementCentroids;

% Assign mu0 (must go above: needed by computeGammaMatrices)
fmdl.mu0 = mu0;

% % Compute and assign dGamma
% for m = 1:length(fmdl.sensors)
%     [dGammaX,dGammaY,dGammaZ] = ...
%         compute_gamma_matrices_derivatives(fmdl,m);
%     dGammaXcell{m} = dGammaX;
%     dGammaYcell{m} = dGammaY;
%     dGammaZcell{m} = dGammaZ;
% end
% 
% fmdl.dGammaXcell = dGammaXcell;
% fmdl.dGammaYcell = dGammaYcell;
% fmdl.dGammaZcell = dGammaZcell;
end