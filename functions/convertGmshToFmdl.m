function fmdl = convertGmshToFmdl(mesh)
% Convert Gmsh mesh from readGmshMesh_v4 to EIDORS fmdl structure

% Ensure EIDORS is in the path
if ~exist('mk_fmdl_from_nodes', 'file')
    error('EIDORS not found in path. Please add EIDORS.');
end

disp('Converting GMSH mesh to EIDORS fmdl')

%% Initialize EIDORS model
fmdl = eidors_obj('fwd_model', 'ImportedGmshMesh');

%% Find tetrahedral elements (type 4)
tetraBlocks = mesh.elements.blocks([mesh.elements.blocks.elementType] == 4);
t_idx = 1;
tets = [];
mat_idx = [];
for blk = tetraBlocks'
    conn = cell2mat(blk.connectivity);
    % Map node tags to row indices in nodes array
    [~, connIdx] = ismember(conn, mesh.nodeTags);
    tets = [tets; connIdx];
    mat_idx = [mat_idx; repmat(blk.entityTag, size(conn, 1), 1)];
    t_idx = t_idx + 1;
end

%% Find triangular surface elements (type 2)
triBlocks = mesh.elements.blocks([mesh.elements.blocks.elementType] == 2);
triangles = [];
tri_mat_idx = [];
for blk = triBlocks'
    conn = cell2mat(blk.connectivity);
    % Map node tags to row indices in nodes array
    [~, connIdx] = ismember(conn, mesh.nodeTags);
    triangles = [triangles; connIdx];
    tri_mat_idx = [tri_mat_idx; repmat(blk.entityTag, size(conn, 1), 1)];
end

%% Set conductivity regions (mappings from regions to values can be defined later)
fmdl.mat_idx = mat_idx;

%% Construct fmdl
fmdl.nodes     = mesh.nodes;  % (nNodes x 2) or (nNodes x 3)
fmdl.elems     = tets;  % (nElems x nNodesPerElem)
fmdl.gnd_node  = 1;  % Scalar node index

%% Set boundary if not already set
if ~isfield(fmdl, 'boundary') || isempty(fmdl.boundary)
    fmdl.boundary = find_boundary(fmdl);
end

%% Create physicalNames map from name => tag (use container because some names in mesh.physicalNames are repeated)
if isfield(mesh, 'physicalNames') && ~isempty(mesh.physicalNames)
    nameMap = containers.Map();
    for i = 1:length(mesh.physicalNames)
        name = mesh.physicalNames(i).name;
        tag  = mesh.physicalNames(i).tag;
        nameMap(name) = tag;
    end
    fmdl.physicalNames = nameMap;
else
    warning('No physical names found in mesh.');
    fmdl.physicalNames = containers.Map(); % create empty map
end

%% MISSING! 

% fmdl.electrode = ...;  % Array of structs with fields `nodes` or `faces`, and `z_contact`

%% Add EIDORS default system matrix calculation method
fmdl = mdl_normalize(fmdl, 1);

%% Complete model
% fmdl = fix_model(fmdl);
end
