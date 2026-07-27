function [J,img] = calc_jacobian_mdeit(img,x,lambdatimesdAdp,A,recon_mode,select_sensor_axis,u_struct)

check_input(A,recon_mode);

if isfield(img,'jacobian')
    J = img.jacobian;
    return;
end

switch recon_mode
    case 'eit'
        J = calc_jacobian_eit(img);

    case 'mdeit1'
        if nargin <7
            J = calc_jacobian_1axis(img,A,select_sensor_axis);
        else
            J = calc_jacobian_1axis(img,A,select_sensor_axis,u_struct);
        end

    case 'mdeit3'
        [Jx,Jy,Jz] = calc_jacobian_3axis(img,A);
        J = [Jx;Jy;Jz];
end

% Always recompute, do no store Jacobian ever
% size_in_megabytes = numel(J)*8/1e6;
% 
% if size_in_megabytes>10
%     warning("Jacobian larger than 10Mb. Skipping storing in img struct")
% else
%     img.jacobian = J;
% end

end

function check_input(A,recon_mode)

% Check if A is function handle or matrix
assert(isa(A, 'function_handle') || isstruct(A),'A must be a function handle or struct');

%MISSING (if A is a struct, check if it has the correct fields. A should
%carry a factorization of the A matrix)

% Check expected struct
if isstruct(A)
    assert(isfield(A,'kind'));
    valid_kinds = {'ldl', 'lu'};
    if ~ismember(A.kind, valid_kinds)
        error('invalid kind "%s". Must be''ldl'', or ''lu''.', A.kind);
    end
    
    switch A.kind
        case 'ldl'
            assert(isfield(A,'L'),'LDL kind must have field L')
            assert(isfield(A,'D'),'LDL kind must have field D')
            assert(isfield(A,'P'),'LDL kind must have field p')
            assert(isfield(A,'n'),'LDL kind must have field n')
        case 'lu'
            assert(isfield(A,'L'),'LU kind must have field L')
            assert(isfield(A,'U'),'LU kind must have field U')
            assert(isfield(A,'pv'),'LU kind must have field pv')
            assert(isfield(A,'qv'),'LU kind must have field qv')
            assert(isfield(A,'n'),'LU kind must have field n')
    end
end

% Check if recon_mode is valid
valid_modes = {'mdeit1', 'mdeit3','eit'};
if ~ismember(recon_mode, valid_modes)
    error('invalid recon_mode "%s". Must be''mdeit1'', or ''mdeit3''.', recon_mode);
end

end


%% FUNCTIONS: calc_jacobian for EIT
function J = calc_jacobian_eit(img)
J = calc_jacobian(img);
end

%% FUNCTIONS:
function J = calc_jacobian_1axis(img,A,select_sensor_axis,u_struct)

if nargin == 4
    error('Not using this now');
end

J = calc_jacobian_1axis_ldl(img,select_sensor_axis);

% if isa(A, 'function_handle')
%     if nargin <4
%         J = calc_jacobian_1axis_A_function_handle(img,A,select_sensor_axis);
%     else
%         J = calc_jacobian_1axis_A_function_handle(img,A,select_sensor_axis,u_struct);
%     end
% elseif isstruct(A)
%     if nargin < 4
%         J = calc_jacobian_1axis_A_factorization_struct(img,A,select_sensor_axis);
%     else
%         J = calc_jacobian_1axis_A_factorization_struct(img,A,select_sensor_axis,u_struct);
%     end
% 
% else
%     error('Unexpected')
% end

end




%% FUNCTIONS:
% These versions are optimized with respect to former ones (avoids
% permutation and solves the adjoint systems with \, not pcg)

function J = calc_jacobian_1axis_A_old(img,A,select_sensor_axis)

mu0 = img.fwd_model.mu0;

n_nodes = size(img.fwd_model.nodes,1);
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
    case 2
        Gamma = img.Gamma2;
    case 3
        Gamma = img.Gamma3;
    otherwise
        error('here')
end

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% Solve the adjoint problem for each sensor to get lambda vectors
GammaT = Gamma.';

A_matrix = A(img.elem_data);
reciprocal_cond_estimate = 1 / condest(A_matrix);

% If matrix is ill-conditioned, fall back to pcg method
if reciprocal_cond_estimate < 1e-15 % Numerically singular
    
    lambda = zeros(n_nodes,num_sensors);
    
    % Jacobi preconditioner - matrix free
    d = sqrt(diag(A_matrix));        % vector of diagonal entries

    Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
    Nfun = @(x) x ./ d;              % right preconditioner

    tol = 1e-10;
    num_elements = numel(img.elem_data);
    for m = 1:num_sensors
        [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),tol,num_elements,Mfun,Nfun);
    end

else
    lambda = A_matrix \ (-GammaT);
end

Gx_times_lambda = G.Gx*lambda;
Gy_times_lambda = G.Gy*lambda;
Gz_times_lambda = G.Gz*lambda;

Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

% We want to broadcast arrays into num_sensors*num_stim*num_elems, in that order, so we avoid a permute in dfd
% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [1 1 n_elem]);

GxL = reshape(Gx_times_lambda.', [num_sensors 1 n_elem]); % [: × 1 × numSensors]
GyL = reshape(Gy_times_lambda.', [num_sensors 1 n_elem]);
GzL = reshape(Gz_times_lambda.', [num_sensors 1 n_elem]);

% Expand u-terms to 3D
GxU = reshape(Gx_times_u.', [1 num_stim n_elem]); % [: × numStim × 1]
GyU = reshape(Gy_times_u.', [1 num_stim n_elem]);
GzU = reshape(Gz_times_u.', [1 num_stim n_elem]);

% Compute all dfdx for all sensors+stim
dfdx = elemV .* ( ...
    GxL.*GxU + ...
    GyL.*GyU + ...
    GzL.*GzU );

% Compute all dfdp (also 3D)
Rx_ = reshape(R.Rx, [num_sensors 1 n_elem]);
Ry_ = reshape(R.Ry, [num_sensors 1 n_elem]);
Rz_ = reshape(R.Rz, [num_sensors 1 n_elem]);

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
gx = reshape(g(:,select_sensor_axis,1), [num_sensors 1 1 ]);
gy = reshape(g(:,select_sensor_axis,2), [num_sensors 1 1]);
gz = reshape(g(:,select_sensor_axis,3), [num_sensors 1 1]);

dfdp = mu_factor*(...
    gx.*dCxdp +...
    gy.*dCydp +...
    gz.*dCzdp);

dfd = dfdx + dfdp;   % size: [numStim × numSensors × numElems]

% Now reshape to match J(ids,:)

% collapse first 2 dims → [numSensors*numStim × numElems]
J = reshape(dfd, num_sensors*num_stim, n_elem);

return
end

function J = calc_jacobian_1axis_ldl(img,select_sensor_axis)

mu0 = img.fwd_model.mu0;

n_nodes = size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);
num_electrodes = numel(img.fwd_model.electrode);

% Compute Gamma matrices
img = compute_gamma_matrices(img);

R = img.fwd_model.R;
G = img.fwd_model.G;

switch select_sensor_axis
    case 1
        Gamma = img.Gamma1;
    case 2
        Gamma = img.Gamma2;
    case 3
        Gamma = img.Gamma3;
    otherwise
        error('here')
end

% Factorize lhs system matrix
A_matrix = lhs_eit_full(img);
F = factorise_symmetric(A_matrix);

% Compute EIT forward solution for each current injection pattern
I = zeros(num_electrodes,num_stim);
for j = 1:num_stim
    I(:,j) = img.fwd_model.stimulation(j).stim_pattern;
end

rhs = sparse(n_nodes+num_electrodes,num_stim);
rhs(end-num_electrodes+1:end,:) = I;

u = solve_fact_multiple_rhs(F,rhs);
u = u(1:end-num_electrodes,:);

% Solve the adjoint problem for each sensor to get lambda vectors
GammaT = Gamma.';

rhs = sparse(n_nodes+num_electrodes,num_sensors);
rhs(1:end-num_electrodes,:) = -GammaT;
lambda = solve_fact_multiple_rhs(F,rhs);
lambda = lambda(1:end-num_electrodes,:);

Gx_times_lambda = G.Gx*lambda;
Gy_times_lambda = G.Gy*lambda;
Gz_times_lambda = G.Gz*lambda;

Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

% We want to broadcast arrays into num_sensors*num_stim*num_elems, in that order, so we avoid a permute in dfd
% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [1 1 n_elem]);

GxL = reshape(Gx_times_lambda.', [num_sensors 1 n_elem]); % [: × 1 × numSensors]
GyL = reshape(Gy_times_lambda.', [num_sensors 1 n_elem]);
GzL = reshape(Gz_times_lambda.', [num_sensors 1 n_elem]);

% Expand u-terms to 3D
GxU = reshape(Gx_times_u.', [1 num_stim n_elem]); % [: × numStim × 1]
GyU = reshape(Gy_times_u.', [1 num_stim n_elem]);
GzU = reshape(Gz_times_u.', [1 num_stim n_elem]);

% Compute all dfdx for all sensors+stim
dfdx = elemV .* ( ...
    GxL.*GxU + ...
    GyL.*GyU + ...
    GzL.*GzU );

% Compute all dfdp (also 3D)
Rx_ = reshape(R.Rx, [num_sensors 1 n_elem]);
Ry_ = reshape(R.Ry, [num_sensors 1 n_elem]);
Rz_ = reshape(R.Rz, [num_sensors 1 n_elem]);

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
gx = reshape(g(:,select_sensor_axis,1), [num_sensors 1 1 ]);
gy = reshape(g(:,select_sensor_axis,2), [num_sensors 1 1]);
gz = reshape(g(:,select_sensor_axis,3), [num_sensors 1 1]);

dfdp = mu_factor*(...
    gx.*dCxdp +...
    gy.*dCydp +...
    gz.*dCzdp);

dfd = dfdx + dfdp;   % size: [numStim × numSensors × numElems]

% Now reshape to match J(ids,:)

% collapse first 2 dims → [numSensors*numStim × numElems]
J = reshape(dfd, num_sensors*num_stim, n_elem);

return
end

function J = calc_jacobian_1axis_A_function_handle(img,A,select_sensor_axis,u_struct)

if nargin == 4
    error('Unexpected parameter');
end

mu0 = img.fwd_model.mu0;

n_nodes = size(img.fwd_model.nodes,1);
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
    case 2
        Gamma = img.Gamma2;
    case 3
        Gamma = img.Gamma3;
    otherwise
        error('here')
end

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% Solve the adjoint problem for each sensor to get lambda vectors
GammaT = Gamma.';

A_matrix = A(img.elem_data);

% reciprocal_cond_estimate = 1 / condest(A_matrix); % this causes out of memory error if large number of elements!
% reciprocal_cond_estimate = rcond(A_matrix); % this is supposed to work better
% 
% % If matrix is ill-conditioned, fall back to pcg method
% if reciprocal_cond_estimate < 1e-15 % Numerically singular
% 
%     lambda = zeros(n_nodes,num_sensors);
% 
%     % Jacobi preconditioner - matrix free
%     d = sqrt(diag(A_matrix));        % vector of diagonal entries
% 
%     Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
%     Nfun = @(x) x ./ d;              % right preconditioner
% 
%     tol = 1e-10;
%     num_elements = numel(img.elem_data);
%     for m = 1:num_sensors
%         [lambda(:,m),~,~] = pcg(A_matrix,-GammaT(:,m),tol,num_elements,Mfun,Nfun);
%     end
% 
% else

lambda = A_matrix \ (-GammaT);

Gx_times_lambda = G.Gx*lambda;
Gy_times_lambda = G.Gy*lambda;
Gz_times_lambda = G.Gz*lambda;

Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

% We want to broadcast arrays into num_sensors*num_stim*num_elems, in that order, so we avoid a permute in dfd
% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [1 1 n_elem]);

GxL = reshape(Gx_times_lambda.', [num_sensors 1 n_elem]); % [: × 1 × numSensors]
GyL = reshape(Gy_times_lambda.', [num_sensors 1 n_elem]);
GzL = reshape(Gz_times_lambda.', [num_sensors 1 n_elem]);

% Expand u-terms to 3D
GxU = reshape(Gx_times_u.', [1 num_stim n_elem]); % [: × numStim × 1]
GyU = reshape(Gy_times_u.', [1 num_stim n_elem]);
GzU = reshape(Gz_times_u.', [1 num_stim n_elem]);

% Compute all dfdx for all sensors+stim
dfdx = elemV .* ( ...
    GxL.*GxU + ...
    GyL.*GyU + ...
    GzL.*GzU );

% Compute all dfdp (also 3D)
Rx_ = reshape(R.Rx, [num_sensors 1 n_elem]);
Ry_ = reshape(R.Ry, [num_sensors 1 n_elem]);
Rz_ = reshape(R.Rz, [num_sensors 1 n_elem]);

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
gx = reshape(g(:,select_sensor_axis,1), [num_sensors 1 1 ]);
gy = reshape(g(:,select_sensor_axis,2), [num_sensors 1 1]);
gz = reshape(g(:,select_sensor_axis,3), [num_sensors 1 1]);

dfdp = mu_factor*(...
    gx.*dCxdp +...
    gy.*dCydp +...
    gz.*dCzdp);

dfd = dfdx + dfdp;   % size: [numStim × numSensors × numElems]

% Now reshape to match J(ids,:)

% collapse first 2 dims → [numSensors*numStim × numElems]
J = reshape(dfd, num_sensors*num_stim, n_elem);

return
end


function J = calc_jacobian_1axis_A_factorization_struct(img,A,select_sensor_axis,u_struct)
mu0 = img.fwd_model.mu0;

n_nodes = size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_electrodes = numel(img.fwd_model.electrode);
num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);

% Check factorization size : are we working with the operator A u = b, or
% with A (u,U) = (0,I) ? 
if n_nodes == size(A.D,1)
    flag = 'simplified';
    error('Expected full system matrix A(u,U) = (0,I), not simplified version A u = b')
elseif n_nodes + num_electrodes == size(A.D,1)
    flag = 'full';
else
    error('Unexpected')
end

if nargin < 4
    flag2 = 'no_forward_solution';
elseif nargin==4 
    flag2 = 'forward_solution';
else
    error('Here');
end


% Compute Gamma matrices
img = compute_gamma_matrices(img);

R = img.fwd_model.R;
G = img.fwd_model.G;

switch select_sensor_axis
    case 1
        Gamma = img.Gamma1;
    case 2
        Gamma = img.Gamma2;
    case 3
        Gamma = img.Gamma3;
    otherwise
        error('here')
end

% % Compute EIT forward solution for each current injection pattern 
% u = fwd_solve(img);
% u = u.volt;

switch flag2
    case 'no_forward_solution'
        % The forward solution can be computed with the factorization of A as well
        I = zeros(num_electrodes,num_stim);
        for j = 1:num_stim
            I(:,j) = img.fwd_model.stimulation(j).stim_pattern;
        end

        switch flag
            case 'full'
                rhs = sparse(n_nodes+num_electrodes,num_stim);
                rhs(end-num_electrodes+1:end,:) = I;

                u = solve_fact_multiple_rhs(A,rhs);
                u = u(1:end-num_electrodes,:);

            case 'simplified'
                error('Unimplemented')
        end

        Gx_times_u = G.Gx*u;
        Gy_times_u = G.Gy*u;
        Gz_times_u = G.Gz*u;

    case 'forward_solution'

        Gx_times_u = u_struct.Gx_times_u;
        Gy_times_u = u_struct.Gy_times_u;
        Gz_times_u = u_struct.Gz_times_u;
end


% Solve the adjoint problem for each sensor to get lambda vectors
GammaT = Gamma.';

switch flag
    case 'full'

        rhs = sparse(n_nodes+num_electrodes,num_sensors);
        rhs(1:end-num_electrodes,:) = -GammaT;
        lambda = solve_fact_multiple_rhs(A,rhs);
        lambda = lambda(1:end-num_electrodes,:);

    case 'simplified'
        error('Unimplemented')
end

Gx_times_lambda = G.Gx*lambda;
Gy_times_lambda = G.Gy*lambda;
Gz_times_lambda = G.Gz*lambda;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

% We want to broadcast arrays into num_sensors*num_stim*num_elems, in that order, so we avoid a permute in dfd
% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [1 1 n_elem]);

GxL = reshape(Gx_times_lambda.', [num_sensors 1 n_elem]); % [: × 1 × numSensors]
GyL = reshape(Gy_times_lambda.', [num_sensors 1 n_elem]);
GzL = reshape(Gz_times_lambda.', [num_sensors 1 n_elem]);

% Expand u-terms to 3D
GxU = reshape(Gx_times_u.', [1 num_stim n_elem]); % [: × numStim × 1]
GyU = reshape(Gy_times_u.', [1 num_stim n_elem]);
GzU = reshape(Gz_times_u.', [1 num_stim n_elem]);

% Compute all dfdx for all sensors+stim
dfdx = elemV .* ( ...
    GxL.*GxU + ...
    GyL.*GyU + ...
    GzL.*GzU );

% Compute all dfdp (also 3D)
Rx_ = reshape(R.Rx, [num_sensors 1 n_elem]);
Ry_ = reshape(R.Ry, [num_sensors 1 n_elem]);
Rz_ = reshape(R.Rz, [num_sensors 1 n_elem]);

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
gx = reshape(g(:,select_sensor_axis,1), [num_sensors 1 1 ]);
gy = reshape(g(:,select_sensor_axis,2), [num_sensors 1 1]);
gz = reshape(g(:,select_sensor_axis,3), [num_sensors 1 1]);

dfdp = mu_factor*(...
    gx.*dCxdp +...
    gy.*dCydp +...
    gz.*dCzdp);

dfd = dfdx + dfdp;   % size: [numStim × numSensors × numElems]

% Now reshape to match J(ids,:)

% collapse first 2 dims → [numSensors*numStim × numElems]
J = reshape(dfd, num_sensors*num_stim, n_elem);

return
end

function [Jx,Jy,Jz] = calc_jacobian_3axis(img,A)

mu0 = img.fwd_model.mu0;
n_nodes = size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);

% Compute Gamma matrices
img = compute_gamma_matrices(img);

R = img.fwd_model.R;
G = img.fwd_model.G;

Gamma1 = img.Gamma1;
Gamma2 = img.Gamma2;
Gamma3 = img.Gamma3;

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

Gamma1T = Gamma1.';
Gamma2T = Gamma2.';
Gamma3T = Gamma3.';

A_matrix = A(img.elem_data);
reciprocal_cond_estimate = 1 / condest(A_matrix);

% If matrix is ill-conditioned, fall back to pcg method
if reciprocal_cond_estimate < 1e-15 % Numerically singular
    
    lambdaX = zeros(n_nodes,num_sensors);
    lambdaY = zeros(n_nodes,num_sensors);
    lambdaZ = zeros(n_nodes,num_sensors);
    
    % Jacobi preconditioner - matrix free
    d = sqrt(diag(A_matrix));        % vector of diagonal entries

    Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
    Nfun = @(x) x ./ d;              % right preconditioner

    tol = 1e-10;
    num_elements = numel(img.elem_data);

    % Solve the adjoint problem for each sensor to get lambda vectors
    parfor m = 1:num_sensors
        [lambdaX(:,m),~,~] = pcg(A_matrix,-Gamma1T(:,m),tol,num_elements,Mfun,Nfun);
        [lambdaY(:,m),~,~] = pcg(A_matrix,-Gamma2T(:,m),tol,num_elements,Mfun,Nfun);
        [lambdaZ(:,m),~,~] = pcg(A_matrix,-Gamma3T(:,m),tol,num_elements,Mfun,Nfun);
    end
else
    % Solve the adjoint problem for each sensor to get lambda vectors
    lambdaX = A_matrix \ (-Gamma1T);
    lambdaY = A_matrix \ (-Gamma2T);
    lambdaZ = A_matrix \ (-Gamma3T);
end


Gx_times_lambda_X = G.Gx*lambdaX;
Gy_times_lambda_X = G.Gy*lambdaX;
Gz_times_lambda_X = G.Gz*lambdaX;

Gx_times_lambda_Y = G.Gx*lambdaY;
Gy_times_lambda_Y = G.Gy*lambdaY;
Gz_times_lambda_Y = G.Gz*lambdaY;

Gx_times_lambda_Z = G.Gx*lambdaZ;
Gy_times_lambda_Z = G.Gy*lambdaZ;
Gz_times_lambda_Z = G.Gz*lambdaZ;

Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

mu_factor = mu0/(4*pi);

elemV = img.fwd_model.elem_volume(:);      % [numElems × 1]

% We want to broadcast arrays into num_sensors*num_stim*num_elems, in that order, so we avoid a permute in dfd

% Expand elem_volume to cover stim × sensor
elemV = reshape(elemV, [1 1 n_elem]);

% Expand u-terms to 3D
GxU = reshape(Gx_times_u.', [1 num_stim n_elem]); % [: × numStim × 1]
GyU = reshape(Gy_times_u.', [1 num_stim n_elem]);
GzU = reshape(Gz_times_u.', [1 num_stim n_elem]);

% Compute all dfdp (also 3D)
Rx_ = reshape(R.Rx, [num_sensors 1 n_elem]);
Ry_ = reshape(R.Ry, [num_sensors 1 n_elem]);
Rz_ = reshape(R.Rz, [num_sensors 1 n_elem]);

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

for select_sensor_axis = 1:3

    switch select_sensor_axis
        case 1
            % Expand lambda and R terms to 3D
            GxL = reshape(Gx_times_lambda_X.', [num_sensors 1 n_elem]); 
            GyL = reshape(Gy_times_lambda_X.', [num_sensors 1 n_elem]);
            GzL = reshape(Gz_times_lambda_X.', [num_sensors 1 n_elem]);
        case 2
            % Expand lambda and R terms to 3D
            GxL = reshape(Gx_times_lambda_Y.', [num_sensors 1 n_elem]); 
            GyL = reshape(Gy_times_lambda_Y.', [num_sensors 1 n_elem]);
            GzL = reshape(Gz_times_lambda_Y.', [num_sensors 1 n_elem]);
        case 3
            % Expand lambda and R terms to 3D
            GxL = reshape(Gx_times_lambda_Z.', [num_sensors 1 n_elem]); 
            GyL = reshape(Gy_times_lambda_Z.', [num_sensors 1 n_elem]);
            GzL = reshape(Gz_times_lambda_Z.', [num_sensors 1 n_elem]);
    end

    % Compute all dfdx for all sensors+stim
    dfdx = elemV .* ( ...
        GxL.*GxU + ...
        GyL.*GyU + ...
        GzL.*GzU );



    % g: [num_sensors × 3 × 3]
    gx = reshape(g(:,select_sensor_axis,1), [num_sensors 1 1 ]);
    gy = reshape(g(:,select_sensor_axis,2), [num_sensors 1 1]);
    gz = reshape(g(:,select_sensor_axis,3), [num_sensors 1 1]);

    dfdp = mu_factor*(...
        gx.*dCxdp +...
        gy.*dCydp +...
        gz.*dCzdp);

    dfd = dfdx + dfdp;   % size: [numStim × numSensors × numElems]

    % Now reshape to match J(ids,:)

    % collapse first 2 dims → [numSensors*numStim × numElems]
    J = reshape(dfd, num_sensors*num_stim, n_elem);

    switch select_sensor_axis
        case 1
            Jx = J;
        case 2
            Jy = J;
        case 3
            Jz = J;
    end

end

return
end




%% FUNCTIONS
function F = factorise_symmetric(A)
    F.kind = 'ldl';
    try
        [F.L,F.D,F.P] = ldl(A,'vector'); 
        F.n = size(A,1);
    catch
        error('Couldnt do it')
        % [F.L,F.U,F.pv,F.qv] = lu(A,'vector'); 
        % F.kind='lu'; 
        % F.n   = size(A,1);
    end
end

function X = solve_fact_multiple_rhs(F, rhs)

    switch F.kind

        case 'ldl'
            % Permute RHS (each column independently)
            rp = rhs(F.P, :);

            % LDL solves (all column-wise)
            y  = F.L \ rp;
            z  = F.D \ y;
            w  = F.L' \ z;

            % Allocate full solution matrix
            X = zeros(F.n, size(rhs,2));

            % Unpermute rows
            X(F.P, :) = w;

        case 'lu'
            % Row permutation of RHS
            y = rhs(F.pv, :);

            % Triangular solves
            z = F.L \ y;
            w = F.U \ z;

            % Allocate solution
            X = zeros(F.n, size(rhs,2));

            % Column permutation recovery
            X(F.qv, :) = w;

        otherwise
            error('Unknown factorisation kind.');
    end
end
