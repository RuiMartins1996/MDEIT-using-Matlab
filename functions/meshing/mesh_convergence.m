function [fmdl_converged,u_norms,b_norms,n_elems,n_nodes] = ...
    mesh_convergence(model_parameters,model_folder,minsz,maxsz,max_num_of_meshes,mode)

if nargin <5
    mode = 'multiple_stimulation';
end

options = {'no_magnetometers','no_geometry_matrices'};
% 
assert(strcmp(mode,'single_stimulation') || strcmp(mode,'multiple_stimulation'),...
    'Parameter "mode" must either be "single_stimulation" or "multiple_stimulation"');

if isempty(maxsz)
    max_num_of_meshes = 1;
    maxsz = minsz;
end

if strcmp(mode,'multiple_stimulation')
    % Load fmdls once so we can obtain the number of stimulation patterns 
    %[~,fmdls] = mk_mdeit_model(model_parameters,[],options);
    % 
    % u_norms = -inf*ones(numel(fmdls{1}.stimulation),max_num_of_meshes);
    % b_norms = -inf*ones(numel(fmdls{1}.stimulation),max_num_of_meshes);

    % This is expensive, it's better to build u_norms dinamically
    u_norms = [];
    b_norms = [];
else
    u_norms = -inf*ones(1,max_num_of_meshes);
    b_norms = -inf*ones(1,max_num_of_meshes);
end

n_elems = zeros(1,max_num_of_meshes);
n_nodes = zeros(1,max_num_of_meshes);
slopes = zeros(1,max_num_of_meshes-1);

szs = logspace(log10(minsz),log10(maxsz),max_num_of_meshes);

% DEBUG!!!!!!!!!!!!!!!!!!!!
% szs = [6.6943 7.6525 8.7479 10.0000];
% max_num_of_meshes = 4;

id_converged = max_num_of_meshes;

% Loop through forward models until mesh convergence criterion is satisfied
for n = 1:max_num_of_meshes

    fprintf('Mesh %i of %i\n',n,max_num_of_meshes)

    % Make a foward model with given maxsz and save
    model_parameters.maxsz = szs(max_num_of_meshes-n+1);
    [model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder,options);

    fmdl = fmdls{1};

    if strcmp(mode,'multiple_stimulation')
        % Solve homogeneous forward model for all stimulation patterns
    else
        % Solve homogeneous forward model for the first stimulation pattern
        fmdl.stimulation = fmdl.stimulation(1);
    end

    img = mk_image_mdeit(fmdl,1.0);
    data = fwd_solve(img);
    u = data.volt;

    % % Compute the forward MDEIT forward solution
    % data  = fwd_solve_mdeit(img,data);
    % 
    % Bx = data.Bx;
    % By = data.By;
    % Bz = data.Bz;

    % Compute L2-norm of the solution for all stimulation patterns, and
    % take maximum
    for i = 1:numel(fmdl.stimulation)
        u_norms(i,n) = compute_L2_norm(fmdl, u(:,i));
        % b_norms(i,n) = norm([Bx(:,i);By(:,i);Bz(:,i)],2);
    end

    n_elems(n) = size(fmdl.elems,1);
    n_nodes(n) = size(fmdl.nodes,1);
    
    % Function to compute slope of graph of max(norm(u(:,n))) w.r.t. n_elems 
    f = @(n) (max(u_norms(:,n))-max(u_norms(:,n)))/(n_elems(n)-n_elems(n-1));

    if n>1
        slopes(n-1) = f(n);
    
    if convergence_criterion(slopes(n-1))
        id_converged = n;
        break;
    end
    
    end
end

if ~convergence_criterion()
    warning('Convergence criterion was not satisfied. Returning the finest mesh.');
end

% Make/Load a complete forward model with geometry matrices and magnetometers
model_parameters.maxsz = szs(1);

file_name = model_parameters_to_file_name(model_parameters,model_folder);

% If this model already exists, check if it is complete ( has geometry
% matrices and magnetometers assigned)
if not(isempty(file_name)) && exist(file_name, 'file') == 2

    var = load(file_name);
    fmdl = var.fmdl;

    fields_needed = {'G','R','dGammaXcell','dGammaYcell','dGammaZcell','sensors'};

    try
        if not(all(ismember(fields_needed,fieldnames(fmdl))))
            error('This fwd_model does not have all the required fields');
        end
        % If it passes the above assertion, then output fmdl
        fmdl_converged = fmdl;
    catch
        % Recompute if the fwd_model is not complete
        [model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder,{'recompute'});
        fmdl_converged = fmdls{1};

    end
else
    % Compute if the fwd_model file does not exist
    [model_parameters,fmdls] = mk_mdeit_model(model_parameters,model_folder);
    fmdl_converged = fmdls{1};
end

% TODO:
% % Remove unwanted saved forward models
% id_converged = max_num_of_meshes;
% 
% szs_to_remove = szs;
% szs_to_remove(id_converged) = [];
% 
% model_parameters_array = repmat(model_parameters,1,length(szs_to_remove));
% for i = 1:numel(model_parameters_array)
%    model_parameters_array(i).maxsz =  szs_to_remove(i);
% end
% 
% remove_unwanted_saved_models(model_parameters_array,model_folder);

end




% Convergence criterion defined here
function bool = convergence_criterion(slope)
% if slope < 1e-5
%     bool = true;
% else
%     bool = false;
% end
bool = false;
end

function eta = compute_residual_error(fmdl, u, sigma)
% Compute residual-based error indicators for EIDORS tetrahedral mesh

nodes = fmdl.nodes;
elems = fmdl.elems;
nElem = size(elems,1);
eta = zeros(nElem,1);

% Compute element gradients (constant for P1)

[Gx,Gy,Gz,vol,~] = computeGradientMatrix(fmdl);

dux = Gx*u;
duy = Gy*u;
duz = Gz*u;

grad_u = [dux,duy,duz];

% --- Compute flux-jump terms ---
% Build a map of shared faces
faces = sort([elems(:,[1 2 3]);
    elems(:,[1 2 4]);
    elems(:,[1 3 4]);
    elems(:,[2 3 4])], 2);
[uniqueFaces, ~, ic] = unique(faces, 'rows', 'stable');

% Find neighbors
faceCount = accumarray(ic,1);
interiorFaces = uniqueFaces(faceCount==2,:);
boundaryFaces = uniqueFaces(faceCount==1,:);

for f = 1:size(interiorFaces,1)
    % Find adjacent elements
    adjElems = find(ismember(faces, interiorFaces(f,:), 'rows'));
    K1 = ceil(adjElems(1)/4);
    K2 = ceil(adjElems(2)/4);

    % Compute unit normal from K1
    v1 = nodes(interiorFaces(f,2),:) - nodes(interiorFaces(f,1),:);
    v2 = nodes(interiorFaces(f,3),:) - nodes(interiorFaces(f,1),:);
    n = cross(v1, v2); n = n / norm(n);

    % Flux jump
    J = (sigma(K1)*grad_u(K1,:) - sigma(K2)*grad_u(K2,:)) * n';
    hF = sqrt(vol(K1)^(1/3) * vol(K2)^(1/3));
    areaF = norm(cross(v1,v2))/2;

    % Add contribution to both elements
    eta(K1) = eta(K1) + 0.5*hF*(J^2)*areaF;
    eta(K2) = eta(K2) + 0.5*hF*(J^2)*areaF;
end

eta = sqrt(eta);
end

function u_norm = compute_L2_norm(fmdl, u)

elems = fmdl.elems;
nodes = fmdl.nodes;
n_elem = size(elems, 1);

dim = size(fmdl.nodes,2);

I = 0;

if dim ==3
    % local mass matrix for tetrahedron
    M = @(V) (V/20) * (2*eye(4) + ones(4) - eye(4)); 
elseif dim==2
    % local mass matrix for linear triangle
    M = @(A) (A/12) * [2 1 1; 1 2 1; 1 1 2];
else
    error('Unexpected branch');
end

for e = 1:n_elem

    idx = elems(e, :);
    coords = nodes(idx,:);
    f = u(idx);

    if dim == 3
        v1 = coords(1,:); v2 = coords(2,:); v3 = coords(3,:); v4 = coords(4,:);
        V = abs(det([v2-v1; v3-v1; v4-v1]'))/6;
    else 
        % Compute area using the determinant formula
        v1 = coords(1,:); v2 = coords(2,:); v3 = coords(3,:);
        V = abs(det([v2 - v1; v3 - v1])) / 2;
    end
    
    I = I + f' * M(V) * f;
end

u_norm = sqrt(I);

end

