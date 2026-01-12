function fmdl = placeElectrodes(fmdl,electrodesEEGStruct)

radius = 10; % in mm or same units as mesh
electrode_centers = electrodesEEGStruct.electrodes.pos;

num_elec = size(electrode_centers, 1);

fmdl.electrode = [];
for i = 1:size(electrode_centers,1)
    fmdl.electrode = [fmdl.electrode struct('nodes',[],'z_contact',[])];
end

% Get boundary nodes:
boundaryNodesIds = zeros(3*size(fmdl.boundary,1),1);

for i = 1:size(fmdl.boundary,1)
    boundaryNodesIds(3*(i-1)+1:3*i) = fmdl.boundary(i,:);
end

boundaryNodesIds = unique(boundaryNodesIds);

bdX = fmdl.nodes(boundaryNodesIds,1);
bdY = fmdl.nodes(boundaryNodesIds,2);
bdZ = fmdl.nodes(boundaryNodesIds,3);

for i = 1:num_elec
    center = electrode_centers(i, :);
    
    
    % Step 2: Compute distances to other nodes
    dists = sqrt(sum(([bdX,bdY,bdZ]-center).^2,2));
    
    % Step 3: Find triangles within radius
    elecNodeIds= boundaryNodesIds(dists < radius, :);
    
    % X = fmdl.nodes(elecNodeIds,1);
    % Y = fmdl.nodes(elecNodeIds,2);
    % Z = fmdl.nodes(elecNodeIds,3);
    % 
    % S = 5*ones(length(elecNodeIds),1);
    % scatter3(X,Y,Z,S,'k');

    fmdl.electrode(i).nodes = elecNodeIds;
    fmdl.electrode(i).z_contact = 0.01;
end

end

