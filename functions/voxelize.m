function out = voxelize(img,dx)

assert(isfield(img,'fwd_model'));
assert(isfield(img.fwd_model,'nodes'));
assert(isfield(img.fwd_model,'elems'));

% Extract nodes from EIDORS img object
nodes = img.fwd_model.nodes;     % Nnodes × 3

dim = size(nodes,2);

switch dim
    case 2
        out = voxelize_2d(img,dx);
        return;
    case 3
        out = voxelize_3d(img,dx);
        return;
end

end

%% Functions
function out = voxelize_3d(img,dx)

nodes = img.fwd_model.nodes;     % Nnodes × 3

% Compute bounding box
xmin = min(nodes(:,1)); xmax = max(nodes(:,1));
ymin = min(nodes(:,2)); ymax = max(nodes(:,2));
zmin = min(nodes(:,3)); zmax = max(nodes(:,3));

% --- Compute actual span
Lx = xmax - xmin;
Ly = ymax - ymin;
Lz = zmax - zmin;

% --- Number of voxels (round up to fully cover)
Nx = ceil(Lx / dx);
Ny = ceil(Ly / dx);
Nz = ceil(Lz / dx);

% --- Compute adjusted bounding box for symmetric coverage
extra_x = (Nx*dx - Lx)/2;
extra_y = (Ny*dx - Ly)/2;
extra_z = (Nz*dx - Lz)/2;

xmin_adj = xmin - extra_x;
xmax_adj = xmax + extra_x;
ymin_adj = ymin - extra_y;
ymax_adj = ymax + extra_y;
zmin_adj = zmin - extra_z;
zmax_adj = zmax + extra_z;

% --- Define voxel centers (shifted by dx/2)
xv = linspace(xmin_adj + dx/2, xmax_adj - dx/2, Nx);
yv = linspace(ymin_adj + dx/2, ymax_adj - dx/2, Ny);
zv = linspace(zmin_adj + dx/2, zmax_adj - dx/2, Nz);

% --- Generate full 3D grid
[Xv, Yv, Zv] = meshgrid(xv, yv, zv);

%% Interpolate
fprintf('Interpolating into voxel grid\n');

uVoxel = voxelAverageElem3D(img, Xv, Yv, Zv, 5);

out = struct('Xv',Xv,'Yv',Yv,'Zv',Zv,'uVoxel',uVoxel);
end


%% Functions
function out = voxelize_2d(img,dx)

nodes = img.fwd_model.nodes;     % Nnodes × 2

% Compute bounding box
xmin = min(nodes(:,1)); xmax = max(nodes(:,1));
ymin = min(nodes(:,2)); ymax = max(nodes(:,2));

% --- Compute actual span
Lx = xmax - xmin;
Ly = ymax - ymin;

% --- Number of voxels (round up to fully cover)
Nx = ceil(Lx / dx);
Ny = ceil(Ly / dx);

% --- Compute adjusted bounding box for symmetric coverage
extra_x = (Nx*dx - Lx)/2;
extra_y = (Ny*dx - Ly)/2;

xmin_adj = xmin - extra_x;
xmax_adj = xmax + extra_x;
ymin_adj = ymin - extra_y;
ymax_adj = ymax + extra_y;


% --- Define voxel centers (shifted by dx/2)
xv = linspace(xmin_adj + dx/2, xmax_adj - dx/2, Nx);
yv = linspace(ymin_adj + dx/2, ymax_adj - dx/2, Ny);

% --- Generate full 2D grid
[Xv, Yv] = meshgrid(xv, yv);

%% Interpolate
fprintf('Interpolating into voxel grid\n');

uVoxel = voxelAverageElem2D(img, Xv, Yv, 5);

out = struct('Xv',Xv,'Yv',Yv,'uVoxel',uVoxel);
end


%% Functions:

function U_voxel = voxelAverageElem2D(img, Xv, Yv, nSub)
%VOXELAVERAGEELEM Interpolate element-wise FEM data to voxels via sub-voxel sampling.
%
%   U_voxel = voxelAverageElem(img, Xv, Yv, Zv, insideMask, nSub)
%
% Inputs:
%   img         : EIDORS img object with fields
%                   .fwd_model.nodes   (Nnodes x 3)
%                   .fwd_model.elems   (Nelements x 4)
%                   .elem_data         (Nelements x 1)
%   Xv,Yv,Zv    : voxel grid center coordinates (from meshgrid)
%   insideMask  : logical array same size as Xv (true for voxels inside FEM domain)
%   nSub        : number of subvoxels per axis (default 5)
%
% Output:
%   U_voxel     : interpolated voxel values (same size as Xv)

if nargin < 6
    nSub = 5; % default resolution
end

nodes    = img.fwd_model.nodes;
elements = img.fwd_model.elems;
u_elem   = img.elem_data(:);

voxel_coords = [Xv(:), Yv(:)];
ti = tsearchn(nodes, elements, voxel_coords);
insideMask = reshape(~isnan(ti), size(Xv));

[Nx, Ny] = size(Xv);
U_voxel = nan(size(Xv));  % preallocate with NaN

% --- Voxel dimensions
dx = Xv(1,2,1) - Xv(1,1,1);
dy = Yv(2,1,1) - Yv(1,1,1);

% --- Precompute normalized sub-voxel offsets in [-0.5, 0.5]^3
[subx, suby] = ndgrid(linspace(-0.5, 0.5, nSub));
subOffsets = [subx(:), suby(:)];

% --- Find indices of voxels inside domain
[idxInsideI, idxInsideJ] = ind2sub(size(insideMask), find(insideMask));
nInside = numel(idxInsideI);

fprintf('Interpolating %d interior voxels using %dx%d subvoxels...\n', ...
        nInside, nSub, nSub);

for n = 1:nInside
    i = idxInsideI(n);
    j = idxInsideJ(n);

    % Voxel center
    xc = Xv(i,j);
    yc = Yv(i,j);

    % Generate sub-voxel coordinates (in real space)
    subs = [xc, yc] + subOffsets .* [dx, dy];

    % Determine which tetrahedra sub-voxels lie in
    triIndex = tsearchn(nodes, elements, subs);

    % Compute average element value among subvoxels inside
    insideSub = ~isnan(triIndex);
    if any(insideSub)
        U_voxel(i,j) = mean(u_elem(triIndex(insideSub)));
    end
end
end




%% Functions
function U_voxel = voxelAverageElem3D(img, Xv, Yv, Zv, nSub)
%VOXELAVERAGEELEM Interpolate element-wise FEM data to voxels via sub-voxel sampling.
%
%   U_voxel = voxelAverageElem(img, Xv, Yv, Zv, insideMask, nSub)
%
% Inputs:
%   img         : EIDORS img object with fields
%                   .fwd_model.nodes   (Nnodes x 3)
%                   .fwd_model.elems   (Nelements x 4)
%                   .elem_data         (Nelements x 1)
%   Xv,Yv,Zv    : voxel grid center coordinates (from meshgrid)
%   insideMask  : logical array same size as Xv (true for voxels inside FEM domain)
%   nSub        : number of subvoxels per axis (default 5)
%
% Output:
%   U_voxel     : interpolated voxel values (same size as Xv)

if nargin < 6
    nSub = 5; % default resolution
end

nodes    = img.fwd_model.nodes;
elements = img.fwd_model.elems;
u_elem   = img.elem_data(:);

voxel_coords = [Xv(:), Yv(:), Zv(:)];
ti = tsearchn(nodes, elements, voxel_coords);
insideMask = reshape(~isnan(ti), size(Xv));

[Nx, Ny, Nz] = size(Xv);
U_voxel = nan(size(Xv));  % preallocate with NaN

% --- Voxel dimensions
dx = Xv(1,2,1) - Xv(1,1,1);
dy = Yv(2,1,1) - Yv(1,1,1);
dz = Zv(1,1,2) - Zv(1,1,1);

% --- Precompute normalized sub-voxel offsets in [-0.5, 0.5]^3
[subx, suby, subz] = ndgrid(linspace(-0.5, 0.5, nSub));
subOffsets = [subx(:), suby(:), subz(:)];

% --- Find indices of voxels inside domain
[idxInsideI, idxInsideJ, idxInsideK] = ind2sub(size(insideMask), find(insideMask));
nInside = numel(idxInsideI);

fprintf('Interpolating %d interior voxels using %dx%dx%d subvoxels...\n', ...
        nInside, nSub, nSub, nSub);

for n = 1:nInside
    i = idxInsideI(n);
    j = idxInsideJ(n);
    k = idxInsideK(n);

    % Voxel center
    xc = Xv(i,j,k);
    yc = Yv(i,j,k);
    zc = Zv(i,j,k);

    % Generate sub-voxel coordinates (in real space)
    subs = [xc, yc, zc] + subOffsets .* [dx, dy, dz];

    % Determine which tetrahedra sub-voxels lie in
    tetIndex = tsearchn(nodes, elements, subs);

    % Compute average element value among subvoxels inside
    insideSub = ~isnan(tetIndex);
    if any(insideSub)
        U_voxel(i,j,k) = mean(u_elem(tetIndex(insideSub)));
    end
end
end

