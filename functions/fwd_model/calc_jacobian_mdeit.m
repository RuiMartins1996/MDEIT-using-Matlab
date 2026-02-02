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

    % v = ones(length(x),1);
    % 
    % w = ones(numel(img.fwd_model.electrode)*numel(img.fwd_model.sensors),1);
    % 
    % z = vec_jacobian_product_v2(img, x, A, w, select_sensor_axis, verbose);

    J = calc_jacobian_1axis_direct_fully_vectorized(img,x,A,select_sensor_axis,verbose);
    
    % Jv1 = J*v;
    % y = jacobian_vec_product_v2(img, x, A, v, select_sensor_axis, verbose);
    % 
    % if norm(Jv1-y,2)>1e-14
    %     error('hEre');
    % end
    % 
    % for t = 1:20
    %     disp(t)
    %     tic
    %     J = calc_jacobian_1axis_direct_fully_vectorized(img,x,A,select_sensor_axis,verbose);
    %     Jv1 = J*v;
    %     times_1(t) = toc;
    % end
    % 
    % for t = 1:20
    %     disp(t)
    %     tic
    %     y = jacobian_vec_product_v2(img, x, A, v, select_sensor_axis, verbose);
    %     times_2(t) = toc;
    % end
    % 
    % fprintf('time = %2.2f +/- %2.2f\n',mean(times_1),std(times_1));
    % fprintf('time = %2.2f +/- %2.2f\n',mean(times_2),std(times_2));
    % 
    % disp('End');

elseif strcmp(recon_mode,'mdeit3')
    J = calc_jacobian_3axis_direct_fully_vectorized(img,x,A,verbose);
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

%% FUNCTIONS:
function J = calc_jacobian_1axis_direct_fully_vectorized(img,x,A,select_sensor_axis,verbose)
img.elem_data = x;
mu0 = img.fwd_model.mu0;

n_nodes =  size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);

% Compute Gamma matrices
img = compute_gamma_matrices(img);

R = img.fwd_model.R;
G = img.fwd_model.G;

switch select_sensor_axis
    case 1
        Gamma = img.Gamma1;
        R1 = R.Rz';
        R2 = R.Ry';
    case 2
        Gamma = img.Gamma2;
        R1 = R.Rx';
        R2 = R.Rz';
    case 3
        Gamma = img.Gamma3;
        R1 = R.Ry';
        R2 = R.Rx';
    otherwise
        error('here')
end

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% Solve the adjoint problem for each sensor to get lambda vectors
lambda = zeros(n_nodes,num_sensors);

A_matrix = A(x);

% Create a DataQueue to receive progress updates
q = parallel.pool.DataQueue;

% Define what happens when a message from a worker arrives
afterEach(q, @(t) update_progress(t));

    function update_progress(t)
        progress.count = progress.count + 1;
        progress.times(progress.count) = t;

        est_time_left = (num_sensors - progress.count) * mean(progress.times(1:progress.count));

        fprintf('\r ETA %s: %.1f (s)', progress.label, est_time_left);
        % fprintf('%.1f (s)\n', );
        if progress.count == progress.total
            fprintf('Done \n');
        end
    end

times = zeros(num_sensors,1);
% Shared counter and total (accessible inside nested function)
progress = struct('count', 0, 'total', num_sensors,'times',times, 'label', 'solving lambda systems');

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

parfor m = 1:num_sensors
    t_start = tic;
    [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),1e-6,numel(x),Mfun,Nfun);
    % [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),1e-6,numel(x),R,Rt);
    if verbose
        send(q, toc(t_start));  % send elapsed time for this iteration
    end
end

t1 = tic;

Gx_times_lambda = G.Gx*lambda;
Gy_times_lambda = G.Gy*lambda;
Gz_times_lambda = G.Gz*lambda;

Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [n_elem 1 1]);
% Later this will broadcast to [numElems × numStim × numSensors]

% Expand lambda and R terms to 3D
GxL = reshape(Gx_times_lambda, [n_elem 1 num_sensors]); % [: × 1 × numSensors]
GyL = reshape(Gy_times_lambda, [n_elem 1 num_sensors]);
GzL = reshape(Gz_times_lambda, [n_elem 1 num_sensors]);

% Expand u-terms to 3D
GxU = reshape(Gx_times_u, [n_elem num_stim 1]); % [: × numStim × 1]
GyU = reshape(Gy_times_u, [n_elem num_stim 1]);
GzU = reshape(Gz_times_u, [n_elem num_stim 1]);

% Compute all dfdx for all sensors+stim
dfdx = elemV .* ( ...
    GxL.*GxU + ...
    GyL.*GyU + ...
    GzL.*GzU );

% Compute all dfdp (also 3D)
Rx_ = reshape(R.Rx', [n_elem 1 num_sensors]);
Ry_ = reshape(R.Ry', [n_elem 1 num_sensors]);
Rz_ = reshape(R.Rz', [n_elem 1 num_sensors]);

% These are the derivatives with respect to sigma of the C components of the Gamma matrix,
dCxdp = ( -Rz_.*GyU + Ry_.*GzU );
dCydp = ( -Rx_.*GzU + Rz_.*GxU );
dCzdp = ( -Ry_.*GxU + Rx_.*GyU );

% The g matrix does not depend on sigma. 
g = zeros(num_sensors,3,3);
for m = 1:numel(img.fwd_model.sensors)
    g(m,:,:) = [...
        img.fwd_model.sensors(m).axes.axis1;
        img.fwd_model.sensors(m).axes.axis2;
        img.fwd_model.sensors(m).axes.axis3];
end

% g: [num_sensors × 3 × 3]
gx = reshape(g(:,select_sensor_axis,1), [1 1 num_sensors]);
gy = reshape(g(:,select_sensor_axis,2), [1 1 num_sensors]);
gz = reshape(g(:,select_sensor_axis,3), [1 1 num_sensors]);

dfdp = mu_factor*(...
    gx.*dCxdp +...
    gy.*dCydp +...
    gz.*dCzdp);
    
% switch select_sensor_axis
%     case 1
%         dfdp = mu_factor * ( -R1_.*GyU + R2_.*GzU );
% 
%     case 2
%         dfdp = mu_factor * ( -R1_.*GzU + R2_.*GxU );
% 
%     case 3
%         dfdp = mu_factor * ( -R1_.*GxU + R2_.*GyU );
% end

dfd = dfdx + dfdp;   % size: [numElems × numStim × numSensors]

% Now reshape to match J(ids,:)
% permute to [numSensors × numStim × numElems]
dfd = permute(dfd, [3 2 1]);

% collapse first 2 dims → [numSensors*numStim × numElems]
J = reshape(dfd, num_sensors*num_stim, n_elem);

if verbose
    fprintf('\r Done. Took %d (s)\n',toc(t1));
end

return
end

%% FUNCTIONS:
function J = calc_jacobian_3axis_direct_fully_vectorized(img,x,A,verbose)

img.elem_data = x;
mu0 = img.fwd_model.mu0;

n_nodes =  size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);

% Compute Gamma matrices
img = compute_gamma_matrices(img);

R = img.fwd_model.R;
G = img.fwd_model.G;

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% Solve the adjoint problem for each sensor to get lambda vectors
lambda1 = zeros(n_nodes,num_sensors);
lambda2 = zeros(n_nodes,num_sensors);
lambda3 = zeros(n_nodes,num_sensors);

A_matrix = A(x);

% Create a DataQueue to receive progress updates
q = parallel.pool.DataQueue;

% Define what happens when a message from a worker arrives
afterEach(q, @(t) update_progress(t));

    function update_progress(t)
        progress.count = progress.count + 1;
        progress.times(progress.count) = t;

        est_time_left = (num_sensors - progress.count) * mean(progress.times(1:progress.count));

        fprintf('\r ETA %s: %.1f (s)', progress.label, est_time_left);
        % fprintf('%.1f (s)\n', );
        if progress.count == progress.total
            fprintf('Done \n');
        end
    end

times = zeros(num_sensors,1);
% Shared counter and total (accessible inside nested function)
progress = struct('count', 0, 'total', num_sensors,'times',times, 'label', 'solving lambda systems');

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

Gamma1_t = img.Gamma1.';
Gamma2_t = img.Gamma2.';
Gamma3_t = img.Gamma3.';

% Incomplete Cholesky factorization preconditioner seems to be a bit faster
% than Jacobi preconditioner. However, it breaks down when the
% conductivities become negative
% R = ichol(A_matrix);
% Rt = R';

parfor m = 1:num_sensors
    t_start = tic;

    [lambda1(:,m),~,~] = pcg(A_matrix,-Gamma1_t(:,m),1e-6,numel(x),Mfun,Nfun);
    [lambda2(:,m),~,~] = pcg(A_matrix,-Gamma2_t(:,m),1e-6,numel(x),Mfun,Nfun);
    [lambda3(:,m),~,~] = pcg(A_matrix,-Gamma3_t(:,m),1e-6,numel(x),Mfun,Nfun);

    if verbose
        send(q, toc(t_start));  % send elapsed time for this iteration
    end
end

t1 = tic;

Gx_times_lambda1 = G.Gx*lambda1;
Gy_times_lambda1 = G.Gy*lambda1;
Gz_times_lambda1 = G.Gz*lambda1;

Gx_times_lambda2 = G.Gx*lambda2;
Gy_times_lambda2 = G.Gy*lambda2;
Gz_times_lambda2 = G.Gz*lambda2;

Gx_times_lambda3 = G.Gx*lambda3;
Gy_times_lambda3 = G.Gy*lambda3;
Gz_times_lambda3 = G.Gz*lambda3;

Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [n_elem 1 1]);
% Later this will broadcast to [numElems × numStim × numSensors]

% Expand lambda and R terms to 3D
GxL1 = reshape(Gx_times_lambda1, [n_elem 1 num_sensors]); % [: × 1 × numSensors]
GyL1 = reshape(Gy_times_lambda1, [n_elem 1 num_sensors]);
GzL1 = reshape(Gz_times_lambda1, [n_elem 1 num_sensors]);

% Expand lambda and R terms to 3D
GxL2 = reshape(Gx_times_lambda2, [n_elem 1 num_sensors]); % [: × 1 × numSensors]
GyL2 = reshape(Gy_times_lambda2, [n_elem 1 num_sensors]);
GzL2 = reshape(Gz_times_lambda2, [n_elem 1 num_sensors]);

% Expand lambda and R terms to 3D
GxL3 = reshape(Gx_times_lambda3, [n_elem 1 num_sensors]); % [: × 1 × numSensors]
GyL3 = reshape(Gy_times_lambda3, [n_elem 1 num_sensors]);
GzL3 = reshape(Gz_times_lambda3, [n_elem 1 num_sensors]);

% Expand u-terms to 3D
GxU = reshape(Gx_times_u, [n_elem num_stim 1]); % [: × numStim × 1]
GyU = reshape(Gy_times_u, [n_elem num_stim 1]);
GzU = reshape(Gz_times_u, [n_elem num_stim 1]);

% Compute all dfdx for all sensors+stim
dfdx1 = elemV .* ( ...
    GxL1.*GxU + ...
    GyL1.*GyU + ...
    GzL1.*GzU );

dfdx2 = elemV .* ( ...
    GxL2.*GxU + ...
    GyL2.*GyU + ...
    GzL2.*GzU );

dfdx3 = elemV .* ( ...
    GxL3.*GxU + ...
    GyL3.*GyU + ...
    GzL3.*GzU );

% Compute all dfdp (also 3D)
Rx_ = reshape(R.Rx', [n_elem 1 num_sensors]);
Ry_ = reshape(R.Ry', [n_elem 1 num_sensors]);
Rz_ = reshape(R.Rz', [n_elem 1 num_sensors]);

% These are the derivatives with respect to sigma of the C components of the Gamma matrix,
dCxdp = ( -Rz_.*GyU + Ry_.*GzU );
dCydp = ( -Rx_.*GzU + Rz_.*GxU );
dCzdp = ( -Ry_.*GxU + Rx_.*GyU );

% The g matrix does not depend on sigma.
g = zeros(num_sensors,3,3);
for m = 1:numel(img.fwd_model.sensors)
    g(m,:,:) = [...
        img.fwd_model.sensors(m).axes.axis1;
        img.fwd_model.sensors(m).axes.axis2;
        img.fwd_model.sensors(m).axes.axis3];
end

% g: [num_sensors × 3 × 3]
gx1 = reshape(g(:,1,1), [1 1 num_sensors]);
gy1 = reshape(g(:,1,2), [1 1 num_sensors]);
gz1 = reshape(g(:,1,3), [1 1 num_sensors]);

gx2 = reshape(g(:,2,1), [1 1 num_sensors]);
gy2 = reshape(g(:,2,2), [1 1 num_sensors]);
gz2 = reshape(g(:,2,3), [1 1 num_sensors]);

gx3 = reshape(g(:,3,1), [1 1 num_sensors]);
gy3 = reshape(g(:,3,2), [1 1 num_sensors]);
gz3 = reshape(g(:,3,3), [1 1 num_sensors]);

dfdp1 = mu_factor*(gx1.*dCxdp + gy1.*dCydp + gz1.*dCzdp);

dfdp2 = mu_factor*(gx2.*dCydp + gy2.*dCydp + gz2.*dCydp);

dfdp3 = mu_factor*(gx3.*dCzdp + gy3.*dCzdp + gz3.*dCzdp);

df1d = dfdx1 + dfdp1;   % size: [numElems × numStim × numSensors]
df2d = dfdx2 + dfdp2;   % size: [numElems × numStim × numSensors]
df3d = dfdx3 + dfdp3;   % size: [numElems × numStim × numSensors]

% Now reshape to match J(ids,:)
% permute to [numSensors × numStim × numElems]
df1d = permute(df1d, [3 2 1]);
df2d = permute(df2d, [3 2 1]);
df3d = permute(df3d, [3 2 1]);

% collapse first 2 dims → [numSensors*numStim × numElems]
J1 = reshape(df1d, num_sensors*num_stim, n_elem);
J2 = reshape(df2d, num_sensors*num_stim, n_elem);
J3 = reshape(df3d, num_sensors*num_stim, n_elem);

J = [J1;J2;J3];

if verbose
    fprintf('\r Done. Took %d (s)\n',toc(t1));
end

end





%% FUNCTIONS: Jacobian-vector product (matrix-free)
function y = jacobian_vec_product(img, x, A, v, select_sensor_axis, verbose)
% Computes J*v without assembling J
% img               : EIDORS image struct
% x                 : current conductivity vector (numElems x 1)
% A                 : function handle A(x) that returns FEM matrix
% v                 : vector to multiply (numElems x 1)
% select_sensor_axis: 1,2,3
% verbose           : true/false

img.elem_data = x;
mu0 = img.fwd_model.mu0;

numNodes = size(img.fwd_model.nodes,1);
numElems = size(img.fwd_model.elems,1);
numStim = numel(img.fwd_model.stimulation);
numSensors = numel(img.fwd_model.sensors);

assert(size(v,1) == numElems && size(v,2) == 1,'v must be a column vector with number of elements entries');

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
        error('Invalid sensor axis')
end

% Forward solution
u = fwd_solve(img);
u = u.volt;

% Solve adjoint problems for each sensor
lambda = zeros(numNodes, numSensors);
A_matrix = A(x);

d = sqrt(diag(A_matrix));
Mfun = @(x) x./d;
Nfun = @(x) x./d;

GammaT = Gamma.';

parfor m = 1:numSensors
    [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),1e-6,numel(x),Mfun,Nfun);
end

Gx_times_lambda = img.fwd_model.G.Gx*lambda;
Gy_times_lambda = img.fwd_model.G.Gy*lambda;
Gz_times_lambda = img.fwd_model.G.Gz*lambda;

Gx_times_u = img.fwd_model.G.Gx*u;
Gy_times_u = img.fwd_model.G.Gy*u;
Gz_times_u = img.fwd_model.G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

dfdu_v = zeros(numSensors,numStim);
dfdp_v = zeros(numSensors,numStim);

for m = 1:numSensors
    for j = 1:numStim
        % This is all the columns at row (m,j) of the Jacobian matrix
        %(elemV.*Gx_times_lambda(:,m).*Gx_times_u(:,j))

        % The matrix-vector product is row .* vector
        dfdu_v(m,j) = v'*(...
            elemV.*(...
                Gx_times_lambda(:,m).*Gx_times_u(:,j)+...
                Gy_times_lambda(:,m).*Gy_times_u(:,j)+...
                Gz_times_lambda(:,m).*Gz_times_u(:,j)...
                ));

        switch select_sensor_axis
            case 1
                dfdp_v(m,j) = v'*mu_factor*(-R1(:,m).*Gy_times_u(:,j)+ R2(:,m).*Gz_times_u(:,j));
            case 2
                dfdp_v(m,j) = v'*mu_factor*(-R1(:,m).*Gz_times_u(:,j)+ R2(:,m).*Gx_times_u(:,j));
            case 3
                dfdp_v(m,j) = v'*mu_factor*(-R1(:,m).*Gx_times_u(:,j)+ R2(:,m).*Gy_times_u(:,j));
        end

    end
end

% collapse first 2 dims → [numSensors*numStim × numElems]
y = reshape(dfdu_v+dfdp_v, numSensors*numStim, 1);

if verbose
    fprintf('\r Done. Took %d (s)\n',toc(t1));
end

end

function y = jacobian_vec_product_v2(img, x, A, v, select_sensor_axis, verbose)

if verbose, t1 = tic; end

img.elem_data = x;
mu0 = img.fwd_model.mu0;

numNodes   = size(img.fwd_model.nodes,1);
numElems   = size(img.fwd_model.elems,1);
numStim    = numel(img.fwd_model.stimulation);
numSensors = numel(img.fwd_model.sensors);

assert(iscolumn(v) && numel(v)==numElems,'v must be numElems×1');

% --- Gamma matrices
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
        error('Invalid sensor axis')
end

% --- Forward solve
u = fwd_solve(img);
u = u.volt;   % [numNodes × numStim]

% --- Assemble system
A_matrix = A(x);

% Jacobi preconditioner
d = sqrt(diag(A_matrix));
Mfun = @(x) x./d;
Nfun = @(x) x./d;

GammaT = Gamma.';

% --- Adjoint solves
lambda = zeros(numNodes,numSensors);

tol   = 1e-6;
maxit = 200;

parfor m = 1:numSensors
    [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),1e-6,numel(x),Mfun,Nfun);
end

% --- Gradients (elementwise)
G = img.fwd_model.G;

Gx_times_lambda = G.Gx * lambda;   % [numElems × numSensors]
Gy_times_lambda = G.Gy * lambda;
Gz_times_lambda = G.Gz * lambda;

Gxu = G.Gx * u;        % [numElems × numStim]
Gyu = G.Gy * u;
Gzu = G.Gz * u;

% --- Preweight once
elemV = img.fwd_model.elem_volume(:);
w = elemV .* v;        % [numElems × 1]

% --- df/du · v   (vectorized!)
dfdu_v = ...
    (Gx_times_lambda.' .* w.') * Gxu + ...
    (Gy_times_lambda.' .* w.') * Gyu + ...
    (Gz_times_lambda.' .* w.') * Gzu;     % [numSensors × numStim]

% --- df/dp · v
mu_factor = mu0/(4*pi);

switch select_sensor_axis
    case 1
        dfdp_v = mu_factor * ( ...
            (-R1.' .* v.') * Gyu + ...
             ( R2.' .* v.') * Gzu );
    case 2
        dfdp_v = mu_factor * ( ...
            (-R1.' .* v.') * Gzu + ...
             ( R2.' .* v.') * Gxu );
    case 3
        dfdp_v = mu_factor * ( ...
            (-R1.' .* v.') * Gxu + ...
             ( R2.' .* v.') * Gyu );
end

% --- Final vector
y = reshape(dfdu_v + dfdp_v, [], 1);

if verbose
    fprintf('Done. Took %.2f s\n', toc(t1));
end
end

% THIS IS WRONG!!!!!!!!!!!!!
function y = vec_jacobian_product_v2(img, x, A, w, select_sensor_axis, verbose)
% Computes w'*J or J^T*w

img.elem_data = x;
mu0 = img.fwd_model.mu0;

numNodes = size(img.fwd_model.nodes,1);
numElems = size(img.fwd_model.elems,1);
numStim = numel(img.fwd_model.stimulation);
numSensors = numel(img.fwd_model.sensors);

assert(size(w,1) == numStim*numSensors && size(w,2) == 1,'w must be a column vector with numStim*numSensors');

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
        error('Invalid sensor axis')
end

% Forward solution
u = fwd_solve(img);
u = u.volt;

% Solve adjoint problems for each sensor
lambda = zeros(numNodes, numSensors);
A_matrix = A(x);

d = sqrt(diag(A_matrix));
Mfun = @(x) x./d;
Nfun = @(x) x./d;

GammaT = Gamma.';

parfor m = 1:numSensors
    [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),1e-6,numel(x),Mfun,Nfun);
end

w_mat = reshape(w, numSensors,numStim);


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

tensorprod(w_mat,dfdx,[1,2],[3,2])

% collapse first 2 dims → [numSensors*numStim × numElems]
y = reshape(dfdu_v+dfdp_v, numSensors*numStim, 1);

% Now reshape to match J(ids,:)
% permute to [numSensors × numStim × numElems]
temp = permute(dfdx, [3 2 1]);

% collapse first 2 dims → [numSensors*numStim × numElems]
Ju = reshape(temp, numSensors*numStim, numElems);


if verbose
    fprintf('\r Done. Took %d (s)\n',toc(t1));
end
end