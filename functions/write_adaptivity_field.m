function write_adaptivity_field(mdl, filename, targetSizes)
%WRITE_ADAPTIVITY_FIELD Write a .pos file defining target element sizes for Gmsh adaptivity
%
%   write_adaptivity_field(mdl, filename, targetSizes)
%
%   mdl          - EIDORS image or forward model (must contain .nodes and .elems)
%   filename     - output .pos filename, e.g. 'adaptivity.pos'
%   targetSizes  - vector of target element sizes (length = nElems)
%
%   The field will be named "TargetSize" and used by Gmsh's AdaptMesh command.

if nargin < 3
    error('Usage: write_adaptivity_field(mdl, filename, targetSizes)');
end

% Determine input type
if isfield(mdl, 'fwd_model')        % image
    nodes = mdl.fwd_model.nodes;
    elems = mdl.fwd_model.elems;
elseif isfield(mdl, 'nodes') && isfield(mdl, 'elems')  % forward model
    nodes = mdl.nodes;
    elems = mdl.elems;
else
    error('Input must be either an EIDORS image or a forward model containing .nodes and .elems');
end

nElems = size(elems, 1);
assert(length(targetSizes) == nElems, ...
    'targetSizes must have one entry per element.');

fid = fopen(filename, 'w');
if fid < 0
    error('Cannot open file %s for writing.', filename);
end

fprintf(fid, 'View "TargetSize" {\n');

% Compute element centroids
dim = size(nodes, 2);
centroids = zeros(nElems, dim);
for i = 1:nElems
    centroids(i,:) = mean(nodes(elems(i,:),:), 1);
end

% Write each element as a point with the corresponding size
for i = 1:nElems
    if dim == 2
        fprintf(fid, 'SP(%g,%g,0){%g};\n', centroids(i,1), centroids(i,2), targetSizes(i));
    else
        fprintf(fid, 'SP(%g,%g,%g){%g};\n', centroids(i,1), centroids(i,2), centroids(i,3), targetSizes(i));
    end
end

fprintf(fid, '};\n');
fclose(fid);

fprintf('✅ Wrote adaptivity field to "%s" (%d elements)\n', filename, nElems);
end
