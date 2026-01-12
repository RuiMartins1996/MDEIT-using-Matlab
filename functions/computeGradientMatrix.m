%Compute gradient matrix G*u = gradu, element volumes and element centroids
function [Gx,Gy,Gz,V,elementCentroids] = computeGradientMatrix(fmdl)

numElements = size(fmdl.elems,1);
numNodes = size(fmdl.nodes,1);

Gx = sparse(numElements,numNodes);
Gy = sparse(numElements,numNodes);
Gz = sparse(numElements,numNodes);

V = zeros(numElements,1);
elementCentroids = zeros(numElements,3);

for k = 1:numElements
    nodeIds = fmdl.elems(k,:);
    nodes = zeros(3,length(nodeIds));
    for n = 1:length(nodeIds)
        nodes(:,n) = fmdl.nodes(nodeIds(n),:);
    end

    elementCentroids(k,:) = sum(nodes,2)/length(nodeIds);
   
    B = [nodes' ones(4,1)];
    grad1 = B\[1;0;0;0];grad2 = B\[0;1;0;0];grad3 = B\[0;0;1;0];grad4= B\[0;0;0;1];
    grad1 = grad1(1:3);grad2 = grad2(1:3);grad3 = grad3(1:3);grad4 = grad4(1:3);

    V(k) = 1/6*abs(det(B));

    grads = [grad1,grad2,grad3,grad4];

    Gx(k,nodeIds) = Gx(k,nodeIds) + grads(1,:);
    Gy(k,nodeIds) = Gy(k,nodeIds) + grads(2,:);
    Gz(k,nodeIds) = Gz(k,nodeIds) + grads(3,:);
end
end
