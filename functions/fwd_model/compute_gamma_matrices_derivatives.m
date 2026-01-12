% FUNCTION: 
function [dGammaX,dGammaY,dGammaZ] = compute_gamma_matrices_derivatives(mdl,magnet_id)

mu0 = mdl.mu0;

if strcmp(mdl.type,'image')
    fmdl = mdl.fwd_model;
elseif strcmp(mdl.type,'fwd_model')
    fmdl = mdl;
else
    error('Model must be of type image or fwd_model');
end

numElements = length(mdl.elem_volume);
numNodes = size(mdl.G.Gx,2);

dim = size(fmdl.nodes,2);

%At most, the number of non-zeros will be (dim+1) non-zeros per element, because
%each element contributes with entries to only the (dim+1) nodes associated to
%that elemet

idn = zeros(numElements*(dim+1),1);
idk = zeros(numElements*(dim+1),1);
valuesX = zeros(numElements*(dim+1),1);
valuesY = zeros(numElements*(dim+1),1);
valuesZ = zeros(numElements*(dim+1),1);

for k = 1:numElements

    nodeIds = fmdl.elems(k,:);

    dGammaXk = mu0/(4*pi)*(...
        -mdl.R.Rz(magnet_id,k)*mdl.G.Gy(k,nodeIds)+...
        mdl.R.Ry(magnet_id,k)*mdl.G.Gz(k,nodeIds));

    dGammaYk = mu0/(4*pi)*(...
        -mdl.R.Rx(magnet_id,k)*mdl.G.Gz(k,nodeIds)+...
        mdl.R.Rz(magnet_id,k)*mdl.G.Gx(k,nodeIds));

    dGammaZk =  mu0/(4*pi)*(...
        -mdl.R.Ry(magnet_id,k)*mdl.G.Gx(k,nodeIds)+...
        mdl.R.Rx(magnet_id,k)*mdl.G.Gy(k,nodeIds));

    ids = (dim+1)*(k-1)+1:(dim+1)*k;
    idn(ids) = nodeIds;    
    idk(ids) = k*ones((dim+1),1);
    valuesX(ids) = dGammaXk(:);
    valuesY(ids) = dGammaYk(:);
    valuesZ(ids) = dGammaZk(:);

end

dGammaX = sparse(idn,idk,valuesX,numNodes,numElements);
dGammaY = sparse(idn,idk,valuesY,numNodes,numElements);
dGammaZ = sparse(idn,idk,valuesZ,numNodes,numElements);
end


%THIS APPROACH IS WRONG FOR NOW, BUT IM SURE THIS CAN BE MADE FASTER!!!

% function [dGammaX,dGammaY,dGammaZ] = compute_gamma_matrices_derivatives(mdl,magnet_id)
% 
% mu0 = mdl.mu0;
% 
% if strcmp(mdl.type,'image')
%     error('Must be fwd_model')
% elseif strcmp(mdl.type,'fwd_model')
%     fmdl = mdl;
% else
%     error('Model must be of type image or fwd_model');
% end
% 
% numElements = size(fmdl.elems,1);
% numNodes    = size(fmdl.nodes,1);
% dim         = size(fmdl.nodes,2);      % 2 or 3
% NpE         = dim+1;                   % nodes per element
% 
% % ----------------------------------------------------------
% % 1) Build sparse indexing arrays for all elements at once
% % ----------------------------------------------------------
% nodeIds = fmdl.elems;          % K × (dim+1)
% 
% % Vectorize node and elem indices
% idn = nodeIds(:);              % K*(dim+1) × 1
% idk = kron((1:numElements)', ones(NpE,1));  % repeated element index
% 
% % ----------------------------------------------------------
% % 2) Extract G-rows corresponding to element nodes
% % ----------------------------------------------------------
% Gx = mdl.G.Gx;     % K × N
% Gy = mdl.G.Gy;
% Gz = mdl.G.Gz;
% 
% % Gather only the G-values needed: makes NpE×K matrix
% Gx_loc = Gx(sub2ind(size(Gx), idk, idn));
% Gy_loc = Gy(sub2ind(size(Gy), idk, idn));
% Gz_loc = Gz(sub2ind(size(Gz), idk, idn));
% 
% % Reshape to K × (dim+1)
% Gx_loc = reshape(Gx_loc, numElements, NpE);
% Gy_loc = reshape(Gy_loc, numElements, NpE);
% Gz_loc = reshape(Gz_loc, numElements, NpE);
% 
% % ----------------------------------------------------------
% % 3) Evaluate the dGamma formulas for all elements at once
% % ----------------------------------------------------------
% Rx = mdl.R.Rx(magnet_id,:).';   % K×1
% Ry = mdl.R.Ry(magnet_id,:).';
% Rz = mdl.R.Rz(magnet_id,:).';
% 
% coef = mu0/(4*pi);
% 
% dGammaX_vals = coef * ( -Rz .* Gy_loc + Ry .* Gz_loc );
% dGammaY_vals = coef * ( -Rx .* Gz_loc + Rz .* Gx_loc );
% dGammaZ_vals = coef * ( -Ry .* Gx_loc + Rx .* Gy_loc );
% 
% % Flatten values
% valuesX = dGammaX_vals(:);
% valuesY = dGammaY_vals(:);
% valuesZ = dGammaZ_vals(:);
% 
% % ----------------------------------------------------------
% % # 4) Build sparse matrices
% % ----------------------------------------------------------
% dGammaX = sparse(idn, idk, valuesX, numNodes, numElements);
% dGammaY = sparse(idn, idk, valuesY, numNodes, numElements);
% dGammaZ = sparse(idn, idk, valuesZ, numNodes, numElements);
% end