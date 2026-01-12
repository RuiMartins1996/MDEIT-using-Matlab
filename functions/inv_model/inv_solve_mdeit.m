function img = inv_solve_mdeit(imdl, data1, data2,J_pre)
%MY_INV_SOLVE  Custom inverse solver with EIDORS-style interface
%
% Supports absolute/difference reconstruction and multiple modes:
%   'eit'      - standard EIT (uses EIDORS inv_solve)
%   'mdeit1'   - 1-axis MDEIT (uses custom solvers)
%   'mdeit3'   - 3-axis MDEIT (uses custom solvers)
%
% Usage:
%   img = my_inv_solve(imdl, data1)       % absolute
%   img = my_inv_solve(imdl, data1, data2) % difference

fprintf('Starting inverse solve ... \n');
%% --- 1. Determine reconstruction mode --- 
if isfield(imdl, 'recon_mode')
    recon_mode = lower(imdl.recon_mode);
    valid_modes = {'eit', 'mdeit1', 'mdeit3'};
    if ~ismember(recon_mode, valid_modes)
        error('my_inv_solve: invalid recon_mode "%s". Must be ''eit'', ''mdeit1'', or ''mdeit3''.', recon_mode);
    end
else
    recon_mode = 'eit';  % default
end

%% --- 2. If standard EIT, use EIDORS inv_solve and return ---
if strcmp(recon_mode, 'eit')
    if nargin < 3
        img = inv_solve(imdl, data1);
    else
        img = inv_solve(imdl, data1, data2);
    end
    return;
end

%% --- 3. Determine reconstruction type and data ---
% If we're here, then its MDEIT reconstruction

if nargin < 3 || isempty(data2)
    % Absolute reconstruction
    recon_type = 'absolute';
    assert(strcmp(imdl.recon_type, recon_type),'recon_type does not match number of input arguments')
    assert( isvector(data1) && isnumeric(data1),'data1 must be a numerical vector');
    data = data1(:);
else
    % Difference reconstruction
    recon_type = 'difference';
    assert(strcmp(imdl.recon_type,recon_type),'recon_type does not match number of input arguments')

    assert(isvector(data1) && isnumeric(data1) && isvector(data2) && isnumeric(data2),'data1 and data2 must be numerical vectors')

    % Fallback if inputs are raw vectors
    data =  data1(:)-data2(:);
end

if nargin == 4
    assert(~isempty(J_pre),'Precomputed jacobian is empty!');

    n_elem = size(imdl.fwd_model.elems,1);
    n_stim = numel(imdl.fwd_model.stimulation);
    n_sensors = numel(imdl.fwd_model.sensors);
    
    if strcmp(recon_mode,'mdeit1')
        assert(all(size(J_pre) == [n_stim*n_sensors,n_elem]),'Incorrect dimensions for precomputed jacobian');
    elseif strcmp(recon_mode,'mdeit3')
        assert(all(size(J_pre) == [3*n_stim*n_sensors,n_elem]),'Incorrect dimensions for precomputed jacobian');
    else
        error('Should not be here!')
    end
    precomputed_jacobian = true;
else
    precomputed_jacobian = false;
end

%% --- 4. Setup solver parameters ---
[x0, tol, max_iterations, Wsqrt, solver, RtR, lambda,preconditiner_type,jacobian_bkgnd,select_sensor_axis,verbose] = ...
    setup_solver_params(imdl,recon_type);

%% --- 5. Define residual and Jacobian function handles ---
img = mk_image_mdeit(imdl.fwd_model,x0);

if strcmp(recon_type,'absolute')
    res = @(x) calc_residual_mdeit(img, x, data,recon_mode,select_sensor_axis);
elseif strcmp(recon_type,'difference')
    res = @(x) data;
end

lambdatimesdAdp = @(lambda) computeLambdaTimesDaDp(img,lambda);
A = @(sigma) M(img,sigma);

if strcmp(recon_type,'absolute')
    jac = @(x) calc_jacobian_mdeit(img, x,lambdatimesdAdp,A,recon_mode,select_sensor_axis);
elseif strcmp(recon_type,'difference')
    if precomputed_jacobian
        jac = @(x) J_pre;
    else
        J = calc_jacobian_mdeit(img, jacobian_bkgnd,lambdatimesdAdp,A,recon_mode,select_sensor_axis);
        jac = @(x) J;
    end
end

%% --- 6. Call selected custom solver ---
switch lower(solver)
    case 'gn'
        x = gn_solve(res, jac, x0, tol, max_iterations, RtR, lambda,[],verbose);
    case 'lm'
        x = lm_solve(res, jac, x0, tol, max_iterations, RtR, lambda,[],verbose);
    case 'tsvd'
        % Only working for Tikhonov
        if not(strcmp(recon_type,'difference')) && not(strcmp(preconditiner_type,'tikhonov'))
            error('Not implemented yet')
        else
            x = tsvd_solve(res, jac, x0, lambda);
        end
    otherwise
        error('my_inv_solve: unknown solver "%s". Use "gn" or "tsvd".', solver);
end

%& --- 7. Build output image ---

if strcmp(recon_type,'difference')
    img.name = sprintf('Custom difference inverse solution (%s, %s)', upper(solver), recon_mode);
    img.fwd_model = imdl.fwd_model;
    img.elem_data = x-x0;
else
    img.name = sprintf('Custom absolute inverse solution (%s, %s)', upper(solver), recon_mode);
    img.fwd_model = imdl.fwd_model;
    img.elem_data = x;
end

if exist("J")
    img.jacobian = J;
end

end






%% ---------- Helper Functions ----------

% setup_solver_params_custom
function [x0, tol, max_iterations, Wsqrt, solver, RtR, lambda,preconditiner_type,jacobian_bkgnd,select_sensor_axis,verbose] = ...
    setup_solver_params(imdl, recon_type)
%SETUP_SOLVER_PARAMS_CUSTOM Prepare parameters for custom solvers
%
% Returns all solver parameters for MDEIT reconstructions

n_elems = size(imdl.fwd_model.elems, 1);
n_stim = numel(imdl.fwd_model.stimulation);
n_sensors = numel(imdl.fwd_model.sensors);

% Check if valid reconstruction type
allowed_recon_types = {'absolute','difference'};
assert(ismember(recon_type,allowed_recon_types),...
    'Invalid reconstruction type: %s. Allowed types are: %s', ...
    recon_type, strjoin(allowed_recon_types, ', '))

% Decide if status information should be printed
if  isfield(imdl, 'verbose')
    assert(ismember(imdl.verbose,[true,false]),'imdl.verbose must either be true or false');
    verbose = imdl.verbose;
else
    verbose = false;
end

% Get jacobian background conductivity and
% override x0
if isfield(imdl,'jacobian_bkgnd')
    if isstruct(imdl.jacobian_bkgnd)
        jacobian_bkgnd = ones(n_elems,1)*imdl.jacobian_bkgnd.value;
    elseif isvector(imdl.jacobian_bkgnd)
        assert(all(size(imdl.jacobian_bkgnd) == [n_elems,1]),'If imdl.jacobian_bkgnd is a vector, then it must be of size n_elems x 1');
        jacobian_bkgnd = imdl.jacobian_bkgnd;
    else
        error('imdl.jacobian_bkgnd must either be a vector or struct of type struct("value",value)')
    end
else
    error('imdl.jacobian_bkgnd must exist');
end

x0 = jacobian_bkgnd;


% --- tolerance ---
if isfield(imdl, 'tol')
    tol = imdl.tol;
else
    tol = 1e-6;
end

% --- maximum iterations ---
if strcmp(recon_type, 'difference')
    max_iterations = 1;          % difference recon always 1 iteration
else
    if isfield(imdl, 'max_iterations')
        max_iterations = imdl.max_iterations;
    else
        max_iterations = 20;      % default for absolute
    end
end

% --- whitening ---
if isfield(imdl, 'Wsqrt')
    Wsqrt = imdl.Wsqrt;
else
    Wsqrt = 1;                    % default (no whitening)
end

% --- solver type ---

%Check that imdl.solver makes sense
allowed_solvers = {'gn','tsvd','lm'};
if isfield(imdl, 'solver')
    assert(ismember(imdl.solver,allowed_solvers),...
        'Invalid solver: %s. Allowed solvers are: %s', ...
          imdl.solver, strjoin(allowed_solvers, ', '))

        solver = imdl.solver;
else
    solver = 'gn';                % default
end

% --- regularization / preconditioner ---
if isfield(imdl, 'RtR_prior')

    if ~isa(imdl.RtR_prior, 'function_handle')
        error('my_inv_solve: "RtR_prior" must be a function handle.');
    end

    RtR = imdl.RtR_prior;

    % Check if pre-conditioner has correct size
    temp = RtR(ones(n_elems,1),ones(n_stim*n_sensors,n_elems));
    if any(size(temp)~=[n_elems,1])
        error('Preconditioner action on vector does not match the expected dimensions of [n_elems,1]');
    end
    
    % Check if it is Tikhonov
    x_test = rand(n_elems,1);
    if  norm(RtR(x_test,ones(n_stim*n_sensors,n_elems)) - x_test)<1e-12
        preconditiner_type = 'tykhonov';
    else
        preconditiner_type = 'unknown';
    end
else
    RtR = @(xk,jk) speye(numel(xk))*xk;                 % default identity (Tikhonov)
end

% --- hyperparameter ---
if isfield(imdl, 'hyperparameter')
    lambda = imdl.hyperparameter.value;
else %default
    lambda = 1e-7;
end

if isfield(imdl,'select_sensor_axis')
    
    assert(isscalar(imdl.select_sensor_axis),'imdl.select_sensor_axis must be either 1,2 or 3');
    assert(ismember(imdl.select_sensor_axis,[1,2,3]),'imdl.select_sensor_axis must be either 1,2 or 3');

    if ~strcmp(imdl.recon_mode,'mdeit1')
        error('Does not make sense to define imdl.select_sensor_axis for imdl.recon_mode = "mdeit3"')
    else
        select_sensor_axis = imdl.select_sensor_axis;
    end
else %default
    select_sensor_axis = 1;
end

end

function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac-Ae*inv(Ad)*Ae';
end