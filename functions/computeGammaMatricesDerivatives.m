%FUNCTION: 
function [dGammaX,dGammaY,dGammaZ] = computeGammaMatricesDerivatives(mdl,magnetId)

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

%At most, the number of non-zeros will be 4 non-zeros per element, because
%each element contributes with entries to only the four nodes associated to
%that elemet

idn = zeros(numElements*4,1);
idk = zeros(numElements*4,1);
valuesX = zeros(numElements*4,1);
valuesY = zeros(numElements*4,1);
valuesZ = zeros(numElements*4,1);

%This approach provides a performance increase!!!! with 328 elements, this
%function takes 1.05e-2 with this approach and 1.32e-2 with the above
%approach, so around 20% faster. With 557 elements is 1.42e-2 to 1.93e-2,
%about 26%.

for k = 1:numElements
    
    nodeIds = fmdl.elems(k,:);

    RxTimesOmegaNum = mdl.R.Rx(magnetId,k)*mdl.elem_volume(k);
    RyTimesOmegaNum = mdl.R.Ry(magnetId,k)*mdl.elem_volume(k);
    RzTimesOmegaNum = mdl.R.Rz(magnetId,k)*mdl.elem_volume(k);

    dGammaXk = mu0/(4*pi)*(-RzTimesOmegaNum*mdl.G.Gy(k,nodeIds)+RyTimesOmegaNum*mdl.G.Gz(k,nodeIds));
    dGammaYk = mu0/(4*pi)*(-RxTimesOmegaNum*mdl.G.Gz(k,nodeIds)+RzTimesOmegaNum*mdl.G.Gx(k,nodeIds));
    dGammaZk =  mu0/(4*pi)*(-RyTimesOmegaNum*mdl.G.Gx(k,nodeIds)+RxTimesOmegaNum*mdl.G.Gy(k,nodeIds));

    %THE PROBLEM WITH THIS APPROACH is that it's changing the non-zero
    %structure of the sparse matrices dGammaXYZ at each iteration. Thi will
    %be fastif that structure is pre-allocated, and we can do this because
    %the number of nnz will be simply K*4 and we know their indexes
    %(nodeIds,k)

    % dGammaX = dGammaX + ...
    %     sparse(nodeIds,k*ones(length(nodeIds),1),dGammaXk(:),numNodes,numElements);
    % dGammaY = dGammaY + ...
    %     sparse(nodeIds,k*ones(length(nodeIds),1),dGammaYk(:),numNodes,numElements);
    % dGammaZ = dGammaZ + ...
    %     sparse(nodeIds,k*ones(length(nodeIds),1),dGammaZk(:),numNodes,numElements);
    
    ids = 4*(k-1)+1:4*k;
    idn(ids) = nodeIds;    
    idk(ids) = k*ones(4,1);
    valuesX(ids) = dGammaXk(:);
    valuesY(ids) = dGammaYk(:);
    valuesZ(ids) = dGammaZk(:);

end

dGammaX = sparse(idn,idk,valuesX,numNodes,numElements);
dGammaY = sparse(idn,idk,valuesY,numNodes,numElements);
dGammaZ = sparse(idn,idk,valuesZ,numNodes,numElements);
end

