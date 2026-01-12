function [J,img] = calc_jacobian_mdeit(img,x,lambdatimesdAdp,A,recon_mode,select_sensor_axis,verbose)

if nargin <7
    verbose = false;
else
    assert(verbose==true || verbose == false,'Must be boolean');
end

valid_modes = {'mdeit1', 'mdeit3','eit'};

if ~ismember(recon_mode, valid_modes)
    error('my_inv_solve: invalid recon_mode "%s". Must be''mdeit1'', or ''mdeit3''.', recon_mode);
end

if strcmp(recon_mode,'eit')
    J = calc_jacobian_eit(img,x);
    return;
end

assert(isa(A, 'function_handle'),'A must be a function handle');

if isfield(img,'jacobian')
    J = img.jacobian;
    return;
end

if strcmp(recon_mode,'mdeit1')
    % This approach is faster!
    % J = calc_jacobian_1axis_direct(img,x,lambdatimesdAdp,A,select_sensor_axis,verbose);
    J = calc_jacobian_1axis_direct_fully_vectorized(img,x,A,select_sensor_axis,verbose);

elseif strcmp(recon_mode,'mdeit3')
    J = calc_jacobian_3axis_direct(img,x,lambdatimesdAdp,A,verbose);
    % J = calc_jacobian_3axis(img,x,lambdatimesdAdp,A);
else
    error('This branch should never be reached')
end

size_in_megabytes = numel(J)*8/1e6;

if size_in_megabytes>10
    warning("Jacobian larger than 10Mb. Skipping storing in img struct")
else
    img.jacobian = J;
end

end

%% FUNCTIONS: calc_jacobian for EIT
function J = calc_jacobian_eit(img)
img.elem_data = x;
J = calc_jacobian(img);
end

function J = calc_jacobian_1axis_direct(img,x,lambdatimesdAdp,A,select_sensor_axis,verbose)
img.elem_data = x;
mu0 = img.fwd_model.mu0;

numNodes =  size(img.fwd_model.nodes,1);
numElements = size(img.fwd_model.elems,1);

numStim = numel(img.fwd_model.stimulation);
numSensors = numel(img.fwd_model.sensors);

J = zeros(numStim*numSensors,numElements);

% Compute Gamma matrices
img = compute_gamma_matrices(img);

switch select_sensor_axis
    case 1
        Gamma = img.GammaX;
        R1 = img.fwd_model.R.Rz';
        R2 = img.fwd_model.R.Ry';
    case 2
        Gamma = img.GammaY;
        R1 = img.fwd_model.R.Rx';
        R2 = img.fwd_model.R.Rz';
    case 3
        Gamma = img.GammaZ;
        R1 = img.fwd_model.R.Ry';
        R2 = img.fwd_model.R.Rx';
    otherwise
        error('here')
end

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% Solve the adjoint problem for each sensor to get lambda vectors
lambda = zeros(numNodes,numSensors);

A_matrix = A(x);

% Create a DataQueue to receive progress updates
q = parallel.pool.DataQueue;

% Define what happens when a message from a worker arrives
afterEach(q, @(t) update_progress(t));

    function update_progress(t)
        progress.count = progress.count + 1;
        progress.times(progress.count) = t;

        est_time_left = (numSensors - progress.count) * mean(progress.times(1:progress.count));

        fprintf('\r ETA %s: %.1f (s)', progress.label, est_time_left);
        % fprintf('%.1f (s)\n', );
        if progress.count == progress.total
            fprintf('Done \n');
        end
    end

times = zeros(numSensors,1);
% Shared counter and total (accessible inside nested function)
progress = struct('count', 0, 'total', numSensors,'times',times, 'label', 'solving lambda systems');

parfor m = 1:numSensors
    t_start = tic;
    % Simple Jacobi preconditioner
    L = sqrt(diag(diag(A_matrix)));
    U = L;

    [lambda(:,m),~,~] = pcg(A_matrix,-Gamma(m,:)',1e-6,numel(x),L,U);

    if verbose
        send(q, toc(t_start));  % send elapsed time for this iteration
    end
end

% t1 = tic;
% fprintf('Precomputing lambda(m) x dAdp:\n');
% lambdaTimes = cell(numSensors,1);
% for m = 1:numSensors
%     lambdaTimes{m} = lambdatimesdAdp(lambda(:,m));   % sparse(numNodes, numElems)
% end


% Instead of computing lambdaTimesdAdp and then multiplying with ut, lets
% do everything directly

t1 = tic;

Gx_times_lambda = img.fwd_model.G.Gx*lambda;
Gy_times_lambda = img.fwd_model.G.Gy*lambda;
Gz_times_lambda = img.fwd_model.G.Gz*lambda;

Gx_times_u = img.fwd_model.G.Gx*u;
Gy_times_u = img.fwd_model.G.Gy*u;
Gz_times_u = img.fwd_model.G.Gz*u;

mu_factor = mu0/(4*pi);

for m = 1:numSensors

    ids = m:numSensors:numStim*numSensors;

    dfdx_m = img.fwd_model.elem_volume.*(...
        Gx_times_lambda(:,m).*Gx_times_u + ...
        Gy_times_lambda(:,m).*Gy_times_u + ...
        Gz_times_lambda(:,m).*Gz_times_u );
    
    switch select_sensor_axis
        case 1
            dfdp_m = mu_factor*(...
                -R1(:,m).*Gy_times_u+R2(:,m).*Gz_times_u);
        case 2
            dfdp_m = mu_factor*(...
                -R1(:,m).*Gz_times_u+R2(:,m).*Gx_times_u);
        case 3
            dfdp_m = mu_factor*(...
                -R1(:,m).*Gx_times_u+R2(:,m).*Gy_times_u);
    end

    J(ids,:) = (dfdx_m + dfdp_m)';
end

if verbose
    fprintf('\r Done. Took %d (s)\n',toc(t1));
end

return

end

function J = calc_jacobian_1axis_direct_fully_vectorized(img,x,A,select_sensor_axis,verbose)
img.elem_data = x;
mu0 = img.fwd_model.mu0;

numNodes =  size(img.fwd_model.nodes,1);
numElems = size(img.fwd_model.elems,1);

numStim = numel(img.fwd_model.stimulation);
numSensors = numel(img.fwd_model.sensors);

% Compute Gamma matrices
img = compute_gamma_matrices(img);

switch select_sensor_axis
    case 1
        Gamma = img.GammaX;
        R1 = img.fwd_model.R.Rz';
        R2 = img.fwd_model.R.Ry';
    case 2
        Gamma = img.GammaY;
        R1 = img.fwd_model.R.Rx';
        R2 = img.fwd_model.R.Rz';
    case 3
        Gamma = img.GammaZ;
        R1 = img.fwd_model.R.Ry';
        R2 = img.fwd_model.R.Rx';
    otherwise
        error('here')
end

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% Solve the adjoint problem for each sensor to get lambda vectors
lambda = zeros(numNodes,numSensors);

A_matrix = A(x);

% Create a DataQueue to receive progress updates
q = parallel.pool.DataQueue;

% Define what happens when a message from a worker arrives
afterEach(q, @(t) update_progress(t));

    function update_progress(t)
        progress.count = progress.count + 1;
        progress.times(progress.count) = t;

        est_time_left = (numSensors - progress.count) * mean(progress.times(1:progress.count));

        fprintf('\r ETA %s: %.1f (s)', progress.label, est_time_left);
        % fprintf('%.1f (s)\n', );
        if progress.count == progress.total
            fprintf('Done \n');
        end
    end

times = zeros(numSensors,1);
% Shared counter and total (accessible inside nested function)
progress = struct('count', 0, 'total', numSensors,'times',times, 'label', 'solving lambda systems');

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

GammaT = Gamma.';

% Incomplete Cholesky factorization preconditioner seems to be a bit faster 
% than Jacobi preconditioner. However, it breaks down when the
% conductivities become negative
% R = ichol(A_matrix);  
% Rt = R';

parfor m = 1:numSensors
    t_start = tic;
    [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),1e-6,numel(x),Mfun,Nfun);
    % [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),1e-6,numel(x),R,Rt);
    if verbose
        send(q, toc(t_start));  % send elapsed time for this iteration
    end
end

t1 = tic;

Gx_times_lambda = img.fwd_model.G.Gx*lambda;
Gy_times_lambda = img.fwd_model.G.Gy*lambda;
Gz_times_lambda = img.fwd_model.G.Gz*lambda;

Gx_times_u = img.fwd_model.G.Gx*u;
Gy_times_u = img.fwd_model.G.Gy*u;
Gz_times_u = img.fwd_model.G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [numElems 1 1]);  
% Later this will broadcast to [numElems × numStim × numSensors]

% Expand lambda and R terms to 3D
GxL = reshape(Gx_times_lambda, [numElems 1 numSensors]); % [: × 1 × numSensors]
GyL = reshape(Gy_times_lambda, [numElems 1 numSensors]);
GzL = reshape(Gz_times_lambda, [numElems 1 numSensors]);

R1_ = reshape(R1, [numElems 1 numSensors]);
R2_ = reshape(R2, [numElems 1 numSensors]);

% Expand u-terms to 3D
GxU = reshape(Gx_times_u, [numElems numStim 1]); % [: × numStim × 1]
GyU = reshape(Gy_times_u, [numElems numStim 1]);
GzU = reshape(Gz_times_u, [numElems numStim 1]);

% Compute all dfdx for all sensors+stim
dfdx = elemV .* ( ...
       GxL.*GxU + ...
       GyL.*GyU + ...
       GzL.*GzU );

% Compute all dfdp (also 3D)
switch select_sensor_axis
    case 1
        dfdp = mu_factor * ( -R1_.*GyU + R2_.*GzU );

    case 2
        dfdp = mu_factor * ( -R1_.*GzU + R2_.*GxU );

    case 3
        dfdp = mu_factor * ( -R1_.*GxU + R2_.*GyU );
end

dfd = dfdx + dfdp;   % size: [numElems × numStim × numSensors]

% Now reshape to match J(ids,:)
% permute to [numSensors × numStim × numElems]
dfd = permute(dfd, [3 2 1]);

% collapse first 2 dims → [numSensors*numStim × numElems]
J = reshape(dfd, numSensors*numStim, numElems);

if verbose
    fprintf('\r Done. Took %d (s)\n',toc(t1));
end

return
end
%% FUNCTIONS:  calc_jacobian_1_axis
function J = calc_jacobian_1axis(img,x,lambdatimesdAdp,A,select_sensor_axis)
img.elem_data = x;

numNodes =  size(img.fwd_model.nodes,1);
numElements = size(img.fwd_model.elems,1);

numStim = numel(img.fwd_model.stimulation);
numSensors = numel(img.fwd_model.sensors);

% img.fwd_model

J = zeros(numStim*numSensors,numElements);

% Compute Gamma matrices
img = compute_gamma_matrices(img);

switch select_sensor_axis
    case 1
        Gamma = img.GammaX;
        dGammaCell = img.fwd_model.dGammaXcell;
    case 2
        Gamma = img.GammaY;
        dGammaCell = img.fwd_model.dGammaYcell;
    case 3
        Gamma = img.GammaZ;
        dGammaCell = img.fwd_model.dGammaZcell;
    otherwise
        error('here')
end

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% Solve the adjoint problem for each sensor to get lambda vectors
lambda = zeros(numNodes,numSensors);

A_matrix = A(x);

% Create a DataQueue to receive progress updates
q = parallel.pool.DataQueue;

% Define what happens when a message from a worker arrives
afterEach(q, @(t) update_progress(t));

    function update_progress(t)
        progress.count = progress.count + 1;
        progress.times(progress.count) = t;

        est_time_left = (numSensors - progress.count) * mean(progress.times(1:progress.count));

        fprintf('\r ETA %s: %.1f (s)', progress.label, est_time_left);
        % fprintf('%.1f (s)\n', );
        if progress.count == progress.total
            fprintf('Done \n');
        end
    end

times = zeros(numSensors,1);
% Shared counter and total (accessible inside nested function)
progress = struct('count', 0, 'total', numSensors,'times',times, 'label', 'solving lambda systems');

parfor m = 1:numSensors
    t_start = tic;
    % Simple Jacobi preconditioner
    L = sqrt(diag(diag(A_matrix)));
    U = L;

    [lambda(:,m),~,~] = pcg(A_matrix,-Gamma(m,:)',1e-6,numel(x),L,U);

    send(q, toc(t_start));  % send elapsed time for this iteration
end

% % This ORIGINAL approach takes 37(s). No parfor takes 106(s).
% fprintf('Assembling Jacobian rows:\n');
% times = zeros(numStim,1);
%
% multiple_rows = zeros(numSensors,numElements);
%
% for n = 1:numStim
%     t_start_2 = tic;
%
%     un_transposed = u(:,n)';
%
%     parfor m = 1:numSensors
%         dfdx = un_transposed*lambdatimesdAdp(lambda(:,m));
%
%         dfdp =  un_transposed*dGammaCell{m};
%
%         multiple_rows(m,:) = dfdx + dfdp;
%     end
%     ids = 1+(n-1)*numSensors : numSensors+(n-1)*numSensors;
%     J(ids,:) = multiple_rows;
%     times(n) = toc(t_start_2);
%     est_time_left = (numStim - n) * mean(times(1:n));
%     fprintf('\r ETA %s: %.1f (s)', 'assembling J rows', est_time_left);
% end
% fprintf('Done\n');


% This approach is FASTER than the original approach
% It has the following improvements over the original:
% - Precompute lambdaxdAdp outside parfor/for loop over sensors
% - COMMENT: with this change, this parfor/for loop runs quickly, so parfor
% is slower than for.
% - Vectorize the computation of dfdx and dfdp by directly performing the
% matrix-matrix products u x lambdaTimes{m}
% - COMMENT: Does not make a big difference for the problem size I tested,
% but maybe will make a difference for many sensors and many stims
% - Vectorize computeLambdaTimesDaDp. Before there was a loop over the
% elements, but we can assemble the sparse matrix directly. Placed the code
% for that into the computeLambdaTimesDaDp function in the functions folder
% - COMMENT: This gave the biggest performance increase!

t1 = tic;
fprintf('Precomputing lambda(m) x dAdp:\n');
lambdaTimes = cell(numSensors,1);
for m = 1:numSensors
    lambdaTimes{m} = lambdatimesdAdp(lambda(:,m));   % sparse(numNodes, numElems)
end

fprintf('Assembling Jacobian rows:\n');

ut = u';

for m = 1:numSensors
    %dfdx = ut*lambdaTimes{m};
    %dfdp =  ut*dGammaCell{m};
    ids = m:numSensors:numStim*numSensors;
    J(ids,:) = ut*lambdaTimes{m}+ut*dGammaCell{m};
end

fprintf('\r Done. Took %d (s)\n',toc(t1));





% This approach takes 159 (s)
% fprintf('Computing Jacobian rows (vectorized):\n');
% tic
% for n = 1:numStim
%     un = u(:,n)'; % 1 × numNodes
%
%     % Precompute dfdx and dfdp for all m in one go
%     dfdx_all = zeros(numSensors, numElements);
%     dfdp_all = zeros(numSensors, numElements);
%
%     for m = 1:numSensors
%         dfdx_all(m,:) = un * lambdatimesdAdp(lambda(:,m));
%         dfdp_all(m,:) = un * dGammaCell{m};
%     end
%
%     ids = (1:numSensors) + (n-1)*numSensors;
%     J(ids,:) = dfdx_all + dfdp_all;
% end
% fprintf('Took %i (s)\n',toc);
% fprintf('Done.\n');

% This approach takes 37.35(s), comparable to ORIGINAL approach
% fprintf('Computing Jacobian rows (parallel over stim):\n');
% tic
%
% J_parts = cell(numStim,1);
%
% parfor n = 1:numStim
%     un = u(:,n)';
%     local_rows = zeros(numSensors, numElements);
%
%     for m = 1:numSensors
%         dfdx = un * lambdatimesdAdp(lambda(:,m));
%         dfdp = un * dGammaCell{m};
%         local_rows(m,:) = dfdx + dfdp;
%     end
%
%     J_parts{n} = local_rows; % store each block
% end
%
% % Combine all results after the parfor
% J = vertcat(J_parts{:});
%
% fprintf('Took %.2f (s)\n', toc);
% fprintf('Done\n');

%This approach takes 114(s). Maybe might be better depending on the number
%of elements. INVESTIGATE!
% fprintf('Computing Jacobian on GPU:\n');
% tic
% u_gpu = gpuArray(u);
% J_gpu = gpuArray.zeros(numSensors * numStim, numElements);
%
% for n = 1:numStim
%     un = u_gpu(:,n)'; % 1×N
%     local_rows = gpuArray.zeros(numSensors, numElements);
%     for m = 1:numSensors
%         dfdx = un * lambdatimesdAdp(lambda(:,m)); % must support gpuArray
%         dfdp = un * dGammaCell{m};                 % must support gpuArray
%         local_rows(m,:) = dfdx + dfdp;
%     end
%     ids = (1:numSensors) + (n-1)*numSensors;
%     J_gpu(ids,:) = local_rows;
% end
% J = gather(J_gpu);
% fprintf('Took %i (s)\n',toc);
% fprintf('Done.\n');

return

end

%% FUNCTIONS
function J = calc_jacobian_3axis_direct(img,x,lambdatimesdAdp,A,verbose)

img.elem_data = x;
mu0 = img.fwd_model.mu0;

numNodes =  size(img.fwd_model.nodes,1);
numElements = size(img.fwd_model.elems,1);

numStim = numel(img.fwd_model.stimulation);
numSensors = numel(img.fwd_model.sensors);

J = zeros(3*numStim*numSensors,numElements);

% Compute Gamma matrices
img = compute_gamma_matrices(img);

GammaX = img.GammaX;
GammaY = img.GammaY;
GammaZ = img.GammaZ;

Rxt = img.fwd_model.R.Rx';
Ryt = img.fwd_model.R.Ry';
Rzt = img.fwd_model.R.Rz';

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% Solve the adjoint problem for each sensor to get lambda vectors
lambdaX = zeros(numNodes,numSensors);
lambdaY = zeros(numNodes,numSensors);
lambdaZ = zeros(numNodes,numSensors);

A_matrix = A(x);

% Create a DataQueue to receive progress updates
q = parallel.pool.DataQueue;

% Define what happens when a message from a worker arrives
afterEach(q, @(t) update_progress(t));

    function update_progress(t)
        progress.count = progress.count + 1;
        progress.times(progress.count) = t;

        est_time_left = (numSensors - progress.count) * mean(progress.times(1:progress.count));

        fprintf('\r ETA %s: %.1f (s)', progress.label, est_time_left);
        % fprintf('%.1f (s)\n', );
        if progress.count == progress.total
            fprintf('Done \n');
        end
    end

times = zeros(numSensors,1);
% Shared counter and total (accessible inside nested function)
progress = struct('count', 0, 'total', numSensors,'times',times, 'label', 'solving lambda systems');

parfor m = 1:numSensors

    t_start = tic;

    % Simple Jacobi preconditioner
    L = sqrt(diag(diag(A_matrix)));
    U = L;

    [lambdaX(:,m),~,~] = pcg(A_matrix,-GammaX(m,:)',1e-5,numel(x),L,U);
    [lambdaY(:,m),~,~] = pcg(A_matrix,-GammaY(m,:)',1e-5,numel(x),L,U);
    [lambdaZ(:,m),~,~] = pcg(A_matrix,-GammaZ(m,:)',1e-5,numel(x),L,U);

    if verbose
        send(q, toc(t_start));  % send elapsed time for this iteration
    end
end

t1 = tic;

Gx_times_lambdaX = img.fwd_model.G.Gx*lambdaX;
Gy_times_lambdaX = img.fwd_model.G.Gy*lambdaX;
Gz_times_lambdaX = img.fwd_model.G.Gz*lambdaX;

Gx_times_lambdaY = img.fwd_model.G.Gx*lambdaY;
Gy_times_lambdaY = img.fwd_model.G.Gy*lambdaY;
Gz_times_lambdaY = img.fwd_model.G.Gz*lambdaY;

Gx_times_lambdaZ = img.fwd_model.G.Gx*lambdaZ;
Gy_times_lambdaZ = img.fwd_model.G.Gy*lambdaZ;
Gz_times_lambdaZ = img.fwd_model.G.Gz*lambdaZ;

Gx_times_u = img.fwd_model.G.Gx*u;
Gy_times_u = img.fwd_model.G.Gy*u;
Gz_times_u = img.fwd_model.G.Gz*u;

for m = 1:numSensors

    % Parenthisis are need, because I want to form the index vector and
    % then correct by (d-1)*numSensors*numStim
    ids_x = (m:numSensors:numStim*numSensors) + 0*numSensors*numStim;
    ids_y = (m:numSensors:numStim*numSensors) + 1*numSensors*numStim;
    ids_z = (m:numSensors:numStim*numSensors) + 2*numSensors*numStim;

    dfdx_m = img.fwd_model.elem_volume.*(...
        Gx_times_lambdaX(:,m).*Gx_times_u + ...
        Gy_times_lambdaX(:,m).*Gy_times_u + ...
        Gz_times_lambdaX(:,m).*Gz_times_u );

    dfdpx_m = mu0/(4*pi)*(...
        -Rzt(:,m).*Gy_times_u+Ryt(:,m).*Gz_times_u);
    
    dfdy_m = img.fwd_model.elem_volume.*(...
        Gx_times_lambdaY(:,m).*Gx_times_u + ...
        Gy_times_lambdaY(:,m).*Gy_times_u + ...
        Gz_times_lambdaY(:,m).*Gz_times_u );

    dfdpy_m = mu0/(4*pi)*(...
        -Rxt(:,m).*Gz_times_u+Rzt(:,m).*Gx_times_u);

    dfdz_m = img.fwd_model.elem_volume.*(...
        Gx_times_lambdaZ(:,m).*Gx_times_u + ...
        Gy_times_lambdaZ(:,m).*Gy_times_u + ...
        Gz_times_lambdaZ(:,m).*Gz_times_u );

    dfdpz_m = mu0/(4*pi)*(...
        -Ryt(:,m).*Gx_times_u+Rxt(:,m).*Gy_times_u);

    J(ids_x,:) = (dfdx_m + dfdpx_m)';
    J(ids_y,:) = (dfdy_m + dfdpy_m)';
    J(ids_z,:) = (dfdz_m + dfdpz_m)';
end

if verbose
    fprintf('\r Done. Took %d (s)\n',toc(t1));
end
return
end
%% FUNCTIONS:  calc_jacobian_3_axis
function J = calc_jacobian_3axis(img,x,lambdatimesdAdp,A)

img.elem_data = x;

numNodes =  size(img.fwd_model.nodes,1);
numElements = size(img.fwd_model.elems,1);

numStim = numel(img.fwd_model.stimulation);
numSensors = numel(img.fwd_model.sensors);

J = zeros(3*numStim*numSensors,numElements);

% Compute Gamma matrices
img = compute_gamma_matrices(img);

GammaX = img.GammaX;
GammaY = img.GammaY;
GammaZ = img.GammaZ;

% Get derivative of Gamma matrices
dGammaXcell = img.fwd_model.dGammaXcell;
dGammaYcell = img.fwd_model.dGammaYcell;
dGammaZcell = img.fwd_model.dGammaZcell;

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% Solve the adjoint problem for each sensor to get lambda vectors
lambdaX = zeros(numNodes,numSensors);
lambdaY = zeros(numNodes,numSensors);
lambdaZ = zeros(numNodes,numSensors);

A_matrix = A(x);

% Create a DataQueue to receive progress updates
q = parallel.pool.DataQueue;

% Define what happens when a message from a worker arrives
afterEach(q, @(t) update_progress(t));

    function update_progress(t)
        progress.count = progress.count + 1;
        progress.times(progress.count) = t;

        est_time_left = (numSensors - progress.count) * mean(progress.times(1:progress.count));

        fprintf('\r ETA %s: %.1f (s)', progress.label, est_time_left);
        % fprintf('%.1f (s)\n', );
        if progress.count == progress.total
            fprintf('Done \n');
        end
    end

times = zeros(numSensors,1);
% Shared counter and total (accessible inside nested function)
progress = struct('count', 0, 'total', numSensors,'times',times, 'label', 'solving lambda systems');

parfor m = 1:numSensors

    t_start = tic;

    % Simple Jacobi preconditioner
    L = sqrt(diag(diag(A_matrix)));
    U = L;

    [lambdaX(:,m),~,~] = pcg(A_matrix,-GammaX(m,:)',1e-5,numel(x),L,U);
    [lambdaY(:,m),~,~] = pcg(A_matrix,-GammaY(m,:)',1e-5,numel(x),L,U);
    [lambdaZ(:,m),~,~] = pcg(A_matrix,-GammaZ(m,:)',1e-5,numel(x),L,U);

    send(q, toc(t_start));  % send elapsed time for this iteration

end

% times = zeros(numStim,1);
%
% multiple_rows_x = zeros(numSensors,numElements);
% multiple_rows_y = zeros(numSensors,numElements);
% multiple_rows_z = zeros(numSensors,numElements);
%
% for n = 1:numStim
%     t_start_2 = tic;
%
%     un_transposed = u(:,n)';
%     for m = 1:numSensors
%
%         dfxdu = un_transposed*lambdatimesdAdp(lambdaX(:,m));
%         dfydu = un_transposed*lambdatimesdAdp(lambdaY(:,m));
%         dfzdu = un_transposed*lambdatimesdAdp(lambdaZ(:,m));
%
%         dfxdp =  un_transposed*dGammaXcell{m};
%         dfydp =  un_transposed*dGammaYcell{m};
%         dfzdp =  un_transposed*dGammaZcell{m};
%
%         multiple_rows_x(m,:) = dfxdu + dfxdp;
%         multiple_rows_y(m,:) = dfydu + dfydp;
%         multiple_rows_z(m,:) = dfzdu + dfzdp;
%     end
%
%     idsx = ...
%         1+(n-1)*numSensors:...
%         numSensors+(n-1)*numSensors;
%
%     idsy = ...
%         1+(n-1)*numSensors+1*numSensors*numStim:...
%         numSensors+(n-1)*numSensors+1*numSensors*numStim;
%
%     idsz = ...
%         1+(n-1)*numSensors+2*numSensors*numStim:...
%         numSensors+(n-1)*numSensors+2*numSensors*numStim;
%
%     J(idsx,:) = multiple_rows_x;
%     J(idsy,:) = multiple_rows_y;
%     J(idsz,:) = multiple_rows_z;
%
%     times(n) = toc(t_start_2);
%     est_time_left = (numStim - n) * mean(times(1:n));
%     fprintf('\r ETA %s: %.1f (s)', 'assembling J rows', est_time_left);
% end
% fprintf('Done\n');


t1 = tic;
fprintf('Precomputing lambda(m) x dAdp:\n');
lambdaXTimes = cell(numSensors,1);
lambdaYTimes = cell(numSensors,1);
lambdaZTimes = cell(numSensors,1);

for m = 1:numSensors
    lambdaXTimes{m} = lambdatimesdAdp(lambdaX(:,m));   % sparse(numNodes, numElems)
    lambdaYTimes{m} = lambdatimesdAdp(lambdaY(:,m));   % sparse(numNodes, numElems)
    lambdaZTimes{m} = lambdatimesdAdp(lambdaZ(:,m));   % sparse(numNodes, numElems)
end

fprintf('Assembling Jacobian rows:\n');

ut = u';

for m = 1:numSensors

    dfdx = ut*lambdaXTimes{m};
    dfdy = ut*lambdaYTimes{m};
    dfdz = ut*lambdaZTimes{m};

    dfxdp =  ut*dGammaXcell{m};
    dfydp =  ut*dGammaYcell{m};
    dfzdp =  ut*dGammaZcell{m};

    % Parenthisis are need, because I want to form the index vector and
    % then correct by (d-1)*numSensors*numStim
    ids_x = (m:numSensors:numStim*numSensors) + 0*numSensors*numStim;
    ids_y = (m:numSensors:numStim*numSensors) + 1*numSensors*numStim;
    ids_z = (m:numSensors:numStim*numSensors) + 2*numSensors*numStim;

    J(ids_x,:) = dfdx+dfxdp;
    J(ids_y,:) = dfdy+dfydp;
    J(ids_z,:) = dfdz+dfzdp;
end

fprintf('\r Done. Took %d (s)\n',toc(t1));


return
end

