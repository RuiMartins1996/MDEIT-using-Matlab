function  img_out = optimize_sensor_configuration(...
    imgh,inv_Gamma_prior,inv_Gamma_noise,...
    jac_coor_transf,q_to_x,x_to_q,...
    mode,type,dim,options)

valid_modes = {'a-opt','d-opt'};
assert(ismember(lower(mode),valid_modes),'Valid modes are "a-opt" or "d-opt"');

valid_types = {'mdeit3','mdeit1'};
assert(ismember(lower(type),valid_types),'Valid types are "mdeit3" or "mdeit1"');

% TODO:
% Check if prior and noise covariance have the correct sizes

if nargin<9
    dim = [];
end

if nargin<10
    options = optimoptions('fminunc',...
        'Display','iter','MaxIterations',30,...
        'StepTolerance',1e-5,'OptimalityTolerance',1e-5,...
        'Algorithm','quasi-newton','HessianApproximation','lbfgs',...
        'SpecifyObjectiveGradient',true,'UseParallel',false);

    fprintf( ...
        ['No options provided.\n' ...
        '\t algorithm: fminunc - quasi-newton\n' ...
        '\t hessian: lbfgs\n' ...
        '\t StepTolerance: 1e-5\n' ...
        '\t MaxIterations: 30\n']);

else
    fprintf( ...
        ['Options provided. Using:\n' ...
        '\t algorithm: fminunc - %s \n' ...
        '\t hessian: %s\n' ...
        '\t StepTolerance: %2.2g\n' ...
        '\t MaxIterations: %i\n'],...
        options.Algorithm,options.HessianApproximation,options.StepTolerance,options.MaxIterations);
end

n_sensors = numel(imgh.fwd_model.sensors);
n_stim = numel(imgh.fwd_model.stimulation);
n_elem = size(imgh.fwd_model.elems,1);
n_nodes = size(imgh.fwd_model.nodes,1);

A = @(x) M(imgh,x);

%% Initialize sensor locations
sensor_locations_0 = fetch_sensor_locations(imgh);

%% Map from sensor locations to generalized coordinates and vice-versa

vector_to_sensor_locations_q = @(q) vector_to_sensor_locations_cartesian(q_to_x(q));
sensor_locations_to_vector_q = @(sensor_locations) x_to_q(sensor_locations_to_vector_cartesian(sensor_locations));

%% Define function+gradient function (needed for fminunc) and f_impl and g_impl

f_impl = @(q) f(imgh,q,inv_Gamma_prior,inv_Gamma_noise,A,vector_to_sensor_locations_q,mode,type,dim);
g_impl  = @(q) g(imgh,q,inv_Gamma_prior,inv_Gamma_noise,A,...
    vector_to_sensor_locations_q,jac_coor_transf,mode,type,dim);

function [func,grad] = funcwithgrad(q,f_impl,g_impl)
% Calculate objective f
func = f_impl(q);

if nargout > 1 % gradient required
    grad =  g_impl(q);
end
end

%% Run optimization
fprintf(' %s OED - fminunc\n',mode)

% Initial conditions
q0 = sensor_locations_to_vector_q(sensor_locations_0);

func_grad = @(q) funcwithgrad(q,f_impl,g_impl);

% % Check if gradient is correct with finite differences ( checks out)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g_0 = g_impl(q0);
% 
% grad_fd = zeros(1,numel(q0));
% delta = 1e-3;
% 
% f_0 = f_impl(q0);
% 
% for k = 1:numel(q0)
%     fprintf('Iteration %i/%i\n',k,numel(q0));
%     qnew = q0;
%     qnew(k) = qnew(k) + delta;
% 
%     grad_fd(k) = (f_impl(qnew)-f_0)/delta;
% 
%     if abs(grad_fd(k)-g_0(k))/abs(g_0(k))>0.01
%         warning('Here');
%     end
% end
% disp(grad_fd-g_0)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[q_opt,fval,history] = runfminunc(func_grad,q0,options);

img_out = assign_sensor_locations(imgh,vector_to_sensor_locations_q(q_opt));

%% Report statisticts

%Change in cost function: 
f0 = f_impl(sensor_locations_to_vector_q(sensor_locations_0));
fend = f_impl(q_opt);

fprintf('Change in objective : %2.2f %%\n',(fend-f0)/f0*100)
end



%% FUNCTION: M
% DESCRIPTION: Matrix A to compute the forward solution u, Au = b
function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac-Ae*inv(Ad)*Ae';
end

%% FUNCTION: dphidp_to_dphidq
% DESCRIPTION: Transform gradient of cost function w.r.t. cartesian coordinates to gradient w.r.t. generalized q coordinates
function dphidq = dphidp_to_dphidq(sensor_locations,dphidp,jacobian_coordinate_transformation)

n_sensors = size(sensor_locations,1);

assert(size(jacobian_coordinate_transformation,2) == 3);
assert(size(jacobian_coordinate_transformation,3) == n_sensors);


qmax = size(jacobian_coordinate_transformation,1);

dphidq = zeros(qmax,n_sensors);

for q = 1:qmax
    temp = zeros(1,n_sensors);
    for dim = 1:3
        temp = temp + squeeze(jacobian_coordinate_transformation(q,dim,:)).'.*dphidp(dim,:);
    end
    dphidq(q,:) = temp;
end

end

%% FUNCTIONS: sensor_locations_to_vector_cartesian and vector_to_sensor_locations_cartesian
% DESCRIPTION: Maps the sensor locations array to a column vector of size
% 3*n_sensors with the x-y-z coordinates of the sensors, and vice-versa

% Map sensor locations to vector in cartesian coordinates
function x = sensor_locations_to_vector_cartesian(sensor_locations)

n_sensors = size(sensor_locations,1);

x = zeros(3*n_sensors,1);

x(1:n_sensors) =  sensor_locations(:,1);
x(n_sensors+1:2*n_sensors) = sensor_locations(:,2);
x(2*n_sensors+1:3*n_sensors) = sensor_locations(:,3);

end

% Map vector in cartesian coordinaes to sensor locations
function sensor_locations = vector_to_sensor_locations_cartesian(x)

assert(mod(numel(x),3)==0);
n_sensors = numel(x)/3;

sensor_locations = zeros(n_sensors,3);

sensor_locations(:,1) = x(1:n_sensors);
sensor_locations(:,2) = x(n_sensors+1:2*n_sensors);
sensor_locations(:,3) = x(2*n_sensors+1:3*n_sensors);

end

%% FUNCTION: f 
% DESCRIPTION: Computes the objective function cost for the A-optimality/
% D-optimality objective
function out = f(img,q,inv_Gamma_prior,inv_Gamma_noise,A,vector_to_sensor_locations,opt_mode,mode,dim)

allowed_opt_modes = {'a-opt','d-opt'};
assert(ismember(opt_mode,allowed_opt_modes));

allowed_modes = {'mdeit3','mdeit1'};
assert(ismember(mode,allowed_modes));

if strcmp(mode,'mdeit3') && nargin<9
    dim = 'default';
end

sensor_locations = vector_to_sensor_locations(q);

switch mode
    case 'mdeit1'
        switch opt_mode
            case 'a-opt'
                phi_oed = compute_cost_function_a_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);
            case 'd-opt'
                phi_oed = compute_cost_function_d_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);
        end
    case 'mdeit3'
        switch opt_mode
            case 'a-opt'
                phi_oed = compute_cost_function_a_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);
            case 'd-opt'
                phi_oed = compute_cost_function_d_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);
        end
end

out = phi_oed;
end

%% FUNCTION: g
% DESCRIPTION: Computes the objective function gradient for the A-optimality/
% D-optimality objective w.r.t. the generalized coordinates 
function dphi = g(img,q,inv_Gamma_prior,inv_Gamma_noise,A,...
    vector_to_sensor_locations,jac_coor_transf,opt_mode,mode,dim)

n_sensors = numel(img.fwd_model.sensors);

allowed_opt_modes = {'a-opt','d-opt'};
assert(ismember(opt_mode,allowed_opt_modes));

allowed_modes = {'mdeit3','mdeit1'};
assert(ismember(mode,allowed_modes));

if strcmp(mode,'mdeit3') && nargin<10
    dim = 'default';
end

sensor_locations = vector_to_sensor_locations(q);

switch mode
    case 'mdeit1'
        switch opt_mode
            case 'a-opt'
                dphidp = compute_cost_function_gradient_a_opt_optimized(...
                    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);
            case 'd-opt'
                dphidp = compute_cost_function_gradient_d_opt_optimized(...
                    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim);
        end
    case 'mdeit3'
        switch opt_mode
            case 'a-opt'
                dphidp = compute_cost_function_gradient_a_opt_optimized_3_axis(...
                    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);
            case 'd-opt'
                dphidp = compute_cost_function_gradient_d_opt_optimized_3_axis(...
                    img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A);
        end
end

% Convert cartesian derivatives to other coordinate derivatives (
% generalized)
jacobian_coordinate_transformation = jac_coor_transf(sensor_locations);
dphidq = dphidp_to_dphidq(sensor_locations,dphidp,jacobian_coordinate_transformation);

dphi = reshape(dphidq.', 1, []);

end

%% Function: assign_sensor_locations
% DESCRIPTION: Changes the sensors property of the forward model to assign a new
% location for the sensors
function img = assign_sensor_locations(img,sensor_locations)
assert(numel(img.fwd_model.sensors) == size(sensor_locations,1));
for id = 1: numel(img.fwd_model.sensors)
    img.fwd_model.sensors(id).position = sensor_locations(id,:);
end
end

%% Function: runfminunc(funcgrad,q0,options)
%DESCRIPTION: This works as a wrapper so we can output the internal state
%of fminunc
function [xsol,fval,history] = runfminunc(funcgrad,q0,options)

% Set up shared variables with outfun
history.x = [];
history.fval = [];

options.OutputFcn = @outfun;

    function stop = outfun(x,optimValues,state)
        stop = false;

        switch state
            case 'init'
                history.fval = [history.fval, optimValues.fval];
                history.x = [history.x, x];
            case 'iter'
                history.fval = [history.fval, optimValues.fval];
                history.x = [history.x, x];
            case 'done'
            otherwise
        end
    end

[xsol,fval] = fminunc(funcgrad,q0,options);

end


%% FUNCTIONS: compute_cost_function_a_opt, compute_cost_function_d_opt, compute_cost_function_a_opt_3_axis,compute_cost_function_d_opt_3_axis
% DESCRIPTION: Implementation of cost function computation for different
% mode and type

function cost = compute_cost_function_d_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim)

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
J = calc_jacobian_mdeit_local(img,A,'mdeit1',dim);
% J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

L = chol(H,'lower');
logdetH = 2*sum(log(diag(L)));

cost = -logdetH;
end

function cost = compute_cost_function_a_opt(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A,dim)

n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
J = calc_jacobian_mdeit_local(img,A,'mdeit1',dim);
% J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% CHANGED HERE!
% cost = trace(inv(H));

L = chol(H,'lower');
Hinv = L'\(L\eye(n_elem));
cost = trace(Hinv);
end

function cost = compute_cost_function_d_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A)
% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
J = calc_jacobian_mdeit_local(img,A,'mdeit3');

% [Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A);
% J = [Jx;Jy;Jz];

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

L = chol(H,'lower');
logdetH = 2*sum(log(diag(L)));

cost = -logdetH;
end

function cost = compute_cost_function_a_opt_3_axis(img,sensor_locations,inv_Gamma_prior,inv_Gamma_noise,A)

n_elem = size(img.fwd_model.elems,1);

% Assign sensor locations
img = assign_sensor_locations(img,sensor_locations);

% Compute the jacobian at current sensor locations
J = calc_jacobian_mdeit_local(img,A,'mdeit3');
% [Jx,Jy,Jz] = calc_jacobian_local(img,A);
% J = [Jx;Jy;Jz];

% Define the inverse posterior covariance matrix
H = J.'*inv_Gamma_noise*J+inv_Gamma_prior;

% CHANGED HERE!
% cost = trace(inv(H));

L = chol(H,'lower');
Hinv = L'\(L\eye(n_elem));
cost = trace(Hinv);

end

%% FUNCTIONS: compute_cost_function_gradient_d_opt_optimized,compute_cost_function_gradient_a_opt_optimized,compute_cost_function_gradient_a_opt_optimized_3_axis,compute_cost_function_gradient_d_opt_optimized_3_axis
%DESCRIPTION: Implementation of the cost function gradient computation for
%different mode and type

function dphidp = compute_cost_function_gradient_d_opt_optimized(...
    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

J = calc_jacobian_mdeit_local(img,A,'mdeit1',dim);
% J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);
dJds = compute_dJ_xyz_optimized(img,dim);

H = J.'*Gamma_noise_inv*J + Gamma_prior_inv;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves

% % Solve L'Y = I -> LH^{-1} = Y -> L'W = H^{-1} -> L H^{-2} = W
% Y = L'\eye(n_elem);
% Hinv = L\Y;

Hinv = L'\(L\eye(n_elem));

W = Gamma_noise_inv*J*Hinv;

dphidp = zeros(3,n_sensors);

for p = 1:3
    for m = 1:n_sensors
        ids = m : n_sensors : n_sensors*n_stim;

        dJm = reshape(dJds{p}(m,:,:), [n_stim, n_elem]);
        Wm  = W(ids,:);

        dphidp(p,m) = -2 * sum(Wm(:) .* dJm(:));
    end
end

end

function dphidp = compute_cost_function_gradient_a_opt_optimized(...
    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A,dim)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

J = calc_jacobian_mdeit_local(img,A,'mdeit1',dim);
% J = calc_jacobian_1axis_direct_fully_vectorized_local(img,A,dim);
dJds = compute_dJ_xyz_optimized(img,dim);

H = J.'*Gamma_noise_inv*J + Gamma_prior_inv;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves

% This was wrong
% % Solve L'Y = I -> LH^{-1} = Y -> L'W = H^{-1} -> L H^{-2} = W
% Y = L'\eye(n_elem);
% X = L\Y;
% W = L'\X;
% inv_H2 = L\W;

Hinv = L'\(L\eye(n_elem));
inv_H2 = Hinv*Hinv;

W = Gamma_noise_inv*J*inv_H2;

dphidp = zeros(3,n_sensors);

for p = 1:3
    for m = 1:n_sensors
        ids = m : n_sensors : n_sensors*n_stim;

        dJm = reshape(dJds{p}(m,:,:), [n_stim, n_elem]);
        Wm  = W(ids,:);

        dphidp(p,m) = -2 * sum(Wm(:) .* dJm(:));
    end
end

end

function dphidp = compute_cost_function_gradient_d_opt_optimized_3_axis(...
    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

J = calc_jacobian_mdeit_local(img,A,'mdeit3');
% [Jx,Jy,Jz] = calc_jacobian_3axis_direct_fully_vectorized_local_optimized(img,A);
% J = [Jx;Jy;Jz];

% dJxds2 = compute_dJ_xyz(img,1);
dJ = compute_dJxyz_xyz_optimized(img);

dJxds = cell(1,3);
dJyds = cell(1,3);
dJzds = cell(1,3);

for p = 1:3
    dJxds{p} = dJ{1,p};
    dJyds{p} = dJ{2,p};
    dJzds{p} = dJ{3,p};
end

H = J.'*Gamma_noise_inv*J + Gamma_prior_inv;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves

% % Solve L'Y = I -> LH^{-1} = Y -> L'W = H^{-1} -> L H^{-2} = W
% Y = L'\eye(n_elem);
% Hinv = L\Y;

Hinv = L'\(L\eye(n_elem));

W = Gamma_noise_inv*J*Hinv;

dphidp = zeros(3,n_sensors);

block_size = n_sensors * n_stim;

for p = 1:3
    for m = 1:n_sensors
        % indices inside one block (x, y, or z)
        ids_local = m : n_sensors : block_size;

        % --- X component ---
        ids_x = ids_local;
        dJx_m = reshape(dJxds{p}(m,:,:), [n_stim, n_elem]);
        Wx_m  = W(ids_x,:);
        term_x = sum(Wx_m(:) .* dJx_m(:));

        % --- Y component ---
        ids_y = ids_local + block_size;
        dJy_m = reshape(dJyds{p}(m,:,:), [n_stim, n_elem]);
        Wy_m  = W(ids_y,:);
        term_y = sum(Wy_m(:) .* dJy_m(:));

        % --- Z component ---
        ids_z = ids_local + 2*block_size;
        dJz_m = reshape(dJzds{p}(m,:,:), [n_stim, n_elem]);
        Wz_m  = W(ids_z,:);
        term_z = sum(Wz_m(:) .* dJz_m(:));

        dphidp(p,m) = -2 * (term_x + term_y + term_z);
    end
end
end

function dphidp = compute_cost_function_gradient_a_opt_optimized_3_axis(...
    img,sensor_locations,Gamma_prior_inv,Gamma_noise_inv,A)

n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);
n_elem = size(img.fwd_model.elems,1);

img = assign_sensor_locations(img,sensor_locations);

J = calc_jacobian_mdeit_local(img,A,'mdeit3');

dJ = compute_dJxyz_xyz_optimized(img);

dJxds = cell(1,3);
dJyds = cell(1,3);
dJzds = cell(1,3);

for p = 1:3
    dJxds{p} = dJ{1,p};
    dJyds{p} = dJ{2,p};
    dJzds{p} = dJ{3,p};
end

H = J.'*Gamma_noise_inv*J + Gamma_prior_inv;

L = chol(H,'lower');

% Compute H^{-2} through Cholesky factorization and linear system solves


% % Solve L'Y = I -> LH^{-1} = Y -> L'W = H^{-1} -> L H^{-2} = W
% Y = L'\eye(n_elem);
% X = L\Y;
% W = L'\X;
% inv_H2 = L\W;

Hinv = L'\(L\eye(n_elem));
inv_H2 = Hinv*Hinv;

W = Gamma_noise_inv*J*inv_H2;

dphidp = zeros(3,n_sensors);

block_size = n_sensors * n_stim;

% Precompute index offsets
ids_base = reshape(1:block_size, n_sensors, n_stim);

for p = 1:3

    dJx = squeeze(dJ{1,p});   % [n_sensors x n_stim x n_elem]
    dJy = squeeze(dJ{2,p});
    dJz = squeeze(dJ{3,p});

    for m = 1:n_sensors

        ids_local = ids_base(m,:);

        Wx = W(ids_local, :);
        Wy = W(ids_local + block_size, :);
        Wz = W(ids_local + 2*block_size, :);

        % Use Frobenius inner products (faster than elementwise sum)
        term_x = sum(sum(Wx .* squeeze(dJx(m,:,:))));
        term_y = sum(sum(Wy .* squeeze(dJy(m,:,:))));
        term_z = sum(sum(Wz .* squeeze(dJz(m,:,:))));

        dphidp(p,m) = -2 * (term_x + term_y + term_z);
    end
end


end



%% Compute dJ(xyz) d xyz (NOT SURE IF THIS IS CORRECT STILL!!!!)) 
function [dJ,dJ1,dJ2] = compute_dJ_xyz_optimized(img,dim)
% Computes derivative of the Jacobian w.r.t. sensor positions
% Vectorized over elements and stimulations, loop only over sensors

n_sensors = numel(img.fwd_model.sensors);
n_elem      = size(img.fwd_model.elems,1);
num_stim    = numel(img.fwd_model.stimulation);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

% Element volumes
elemV = img.fwd_model.elem_volume(:);  
elemV = reshape(elemV,[1 1 n_elem]);

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% % G*u for all stim, size: n_elem x num_stim
% GxU = G.Gx*u;
% GyU = G.Gy*u;
% GzU = G.Gz*u;

% G*u for all stim, size: n_elem x num_stim
Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

% Expand u-terms to [1 x n_stim x n_elem]
GxU = reshape(Gx_times_u.', [1 num_stim n_elem]); 
GyU = reshape(Gy_times_u.', [1 num_stim n_elem]);
GzU = reshape(Gz_times_u.', [1 num_stim n_elem]);

% Compute derivative of lambda w.r.t sensor positions (output dR1 and dR2
% to avoid recomputing)

[dlambda,dR1dp,dR2dp] = compute_dlambda_xyz(img,dim); % This is correct


for dim = 1:3
    for p = 1:3
        DLx{p} = reshape((G.Gx * dlambda{p}).',[n_sensors 1 n_elem]);
        DLy{p} = reshape((G.Gy * dlambda{p}).',[n_sensors 1 n_elem]);
        DLz{p} = reshape((G.Gz * dlambda{p}).',[n_sensors 1 n_elem]);
    end
end


% Compute derivatives of R w.r.t sensor positions
switch dim
    case 1
        G1U = GyU;
        G2U = GzU;
    case 2
        G1U = GzU;
        G2U = GxU;
    case 3
        G1U = GxU;
        G2U = GyU;
    otherwise
        error('Dimension not supported.');
end

% Preallocate
dJ = cell(1,3);
dJ1 = zeros(n_sensors,num_stim,n_elem);
dJ2 = zeros(n_sensors,num_stim,n_elem);

%% Loop only over sensors
for p = 1:3

    % --- dJ1 term (fully vectorized over sensors) ---

    tmp = DLx{p} .* GxU + ...
        DLy{p} .* GyU + ...
        DLz{p} .* GzU;   % num_sensors xnum_stim x n_elem

    % Expand element volumes
    dJ1 = tmp .* elemV;

    % --- dJ2 term ---

    term2 = ...
        reshape(dR2dp{p},[n_sensors 1 n_elem]) .* G2U - ...
        reshape(dR1dp{p},[n_sensors 1 n_elem]) .* G1U;

    dJ2 = -mu_factor * term2;

    % Total
    dJ{p} = dJ1 - dJ2;
    
    % for m = 1:n_sensors
    %     %% --- dJ1: contribution from dlambda ---
    %     % dlambda * G matrices, size: n_elem x 1
    %     dlGx = G.Gx*dlambda{p}(:,m);
    %     dlGy = G.Gy*dlambda{p}(:,m);
    %     dlGz = G.Gz*dlambda{p}(:,m);
    % 
    %     % Elementwise multiplication and sum over components
    %     % tmp: n_elem x num_stim
    %     tmp = dlGx .* GxU + dlGy .* GyU + dlGz .* GzU;
    % 
    %     % Multiply by element volumes, permute to [num_stim x n_elem]
    %     dJ1(m,:,:) = tmp.' .* elemV(:).';
    % 
    %     %% --- dJ2: contribution from dR/dp ---
    %     % dRydp, dRzdp: 1 x n_elem
    %     % Multiply by G*u per stimulation, permute to [num_stim x n_elem]
    %     dJ2(m,:,:) = -mu_factor * ( ...
    %         dR2dp{p}(m,:) .* G2U.' - dR1dp{p}(m,:) .* G1U.' ...
    %         );
    % end

    %% Total derivative
    dJ{p} = dJ1 - dJ2;
end

end

function [dJ,dJ1,dJ2] = compute_dJxyz_xyz_optimized(img)
% Computes derivative of the Jacobian w.r.t. sensor positions
% Vectorized over elements and stimulations, loop only over sensors

n_sensors = numel(img.fwd_model.sensors);
n_elem      = size(img.fwd_model.elems,1);
num_stim    = numel(img.fwd_model.stimulation);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

% Element volumes
elemV = img.fwd_model.elem_volume(:);  
elemV = reshape(elemV,[1 1 n_elem]);

% Compute EIT forward solution for each current injection pattern
u = fwd_solve(img);
u = u.volt;

% G*u for all stim, size: n_elem x num_stim
Gx_times_u = G.Gx*u;
Gy_times_u = G.Gy*u;
Gz_times_u = G.Gz*u;

% Expand u-terms to [1 x n_stim x n_elem]
GxU = reshape(Gx_times_u.', [1 num_stim n_elem]); 
GyU = reshape(Gy_times_u.', [1 num_stim n_elem]);
GzU = reshape(Gz_times_u.', [1 num_stim n_elem]);


% Compute derivative of lambda w.r.t sensor positions (output dR1 and dR2
% to avoid recomputing)
[dlambda,dR1dp,dR2dp] = compute_dlambdaxyz_xyz_optimized(img);

% Careful here: dlambda is a 3 x 3 cell of dim x coord. dlambdaXdx =
% dlambda{1,1}
for dim = 1:3
    for p = 1:3
        DLx{dim,p} = reshape((G.Gx * dlambda{dim,p}).',[n_sensors 1 n_elem]);
        DLy{dim,p} = reshape((G.Gy * dlambda{dim,p}).',[n_sensors 1 n_elem]);
        DLz{dim,p} = reshape((G.Gz * dlambda{dim,p}).',[n_sensors 1 n_elem]);
    end
end

% Precompute G1U and G2U
G1U = {GyU, GzU, GxU};
G2U = {GzU, GxU, GyU};

% Preallocate
dJ = cell(3,3);

%% Loop only over sensors
for dim = 1:3
    for p = 1:3
        % --- dJ1 term (fully vectorized over sensors) ---
        

        tmp = DLx{dim,p} .* GxU + ...
            DLy{dim,p} .* GyU + ...
            DLz{dim,p} .* GzU;   % num_sensors xnum_stim x n_elem 

        % THIS WAS WRONG BEFORE ;(. Quite hard to catch the bug!
        % tmp = DLx{p} .* GxU + ...
        %     DLy{p} .* GyU + ...
        %     DLz{p} .* GzU;   % num_sensors xnum_stim x n_elem 

        % Expand element volumes
        dJ1 = tmp .* elemV;

        % --- dJ2 term ---

        term2 = ...
            reshape(dR2dp{dim,p},[n_sensors 1 n_elem]) .* G2U{dim} - ...
            reshape(dR1dp{dim,p},[n_sensors 1 n_elem]) .* G1U{dim};

        dJ2 = -mu_factor * term2;

        % Total
        dJ{dim,p} = dJ1 - dJ2;
    end
end
end

%% Compute dlambda(xyz)d(xyz)
function [dlambda,dR1t,dR2t] = compute_dlambdaxyz_xyz_optimized(img)

num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);

dlambda = cell(3,3);

dlambdaX = zeros(n_nodes,num_sensors);
dlambdaY = zeros(n_nodes,num_sensors);
dlambdaZ = zeros(n_nodes,num_sensors);

dR1t = cell(3,3);
dR2t = cell(3,3);


mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

sigma = img.elem_data;
A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

dR_all{1} = dRmkj_xyz_optimized(img.fwd_model,1);
dR_all{2} = dRmkj_xyz_optimized(img.fwd_model,2);
dR_all{3} = dRmkj_xyz_optimized(img.fwd_model,3);

sigmaT = sigma.';

n_elem = numel(sigma);

for dim = 1:3

    switch dim
        case 1
            dR1 = dR_all{3};
            dR2 = dR_all{2};

            G1 = G.Gy;
            G2 = G.Gz;
        case 2
            dR1  = dR_all{1};
            dR2 = dR_all{3};

            G1 = G.Gz;
            G2 = G.Gx;
        case 3
            dR1 = dR_all{2};
            dR2 = dR_all{1};

            G1 = G.Gx;
            G2 = G.Gy;
    end

    % Multiply elementwise by sigma (diagonal) without creating sparse diagonal matrix
    % rhs = mu_factor*( dR1 * Sigma * G1 -  dR2* Sigma * G2);

    rhs1 = mu_factor*( (dR1{1} .* sigmaT)*G1 - (dR2{1} .* sigmaT)*G2 ); %for p=1
    rhs2 = mu_factor*( (dR1{2} .* sigmaT)*G1 - (dR2{2} .* sigmaT)*G2 ); %for p=2
    rhs3 = mu_factor*( (dR1{3} .* sigmaT)*G1 - (dR2{3} .* sigmaT)*G2 ); %for p=3

    % Solve the adjoint problem for each sensor to get lambda vectors
    parfor m = 1:num_sensors

        [dlambdaX(:,m),~,~] = pcg(A_matrix,rhs1(m,:)',1e-9,n_elem,Mfun,Nfun);
        [dlambdaY(:,m),~,~] = pcg(A_matrix,rhs2(m,:)',1e-9,n_elem,Mfun,Nfun);
        [dlambdaZ(:,m),~,~] = pcg(A_matrix,rhs3(m,:)',1e-9,n_elem,Mfun,Nfun);
    end

    dlambda{dim,1} = dlambdaX;
    dlambda{dim,2} = dlambdaY;
    dlambda{dim,3} = dlambdaZ;

    dR1t{dim,1} = dR1{1};
    dR1t{dim,2} = dR1{2};
    dR1t{dim,3} = dR1{3};

    dR2t{dim,1} = dR2{1};
    dR2t{dim,2} = dR2{2};
    dR2t{dim,3} = dR2{3};
end

end

function [dlambda,dR1,dR2] = compute_dlambda_xyz(img,dim)

num_sensors = numel(img.fwd_model.sensors);
n_nodes = size(img.fwd_model.nodes,1);

dlambda_dx = zeros(n_nodes,num_sensors);
dlambda_dy = zeros(n_nodes,num_sensors);
dlambda_dz = zeros(n_nodes,num_sensors);

mu_factor = img.fwd_model.mu0/(4*pi);

G = img.fwd_model.G;

switch dim
    case 1
        dR1 = dRmkj_xyz(img.fwd_model, 3);
        dR2 = dRmkj_xyz(img.fwd_model, 2);

        G1 = G.Gy;
        G2 = G.Gz;
    case 2
        dR1  = dRmkj_xyz(img.fwd_model, 1);
        dR2 = dRmkj_xyz(img.fwd_model, 3);

        G1 = G.Gz;
        G2 = G.Gx;
    case 3
        dR1 = dRmkj_xyz(img.fwd_model, 2);
        dR2 = dRmkj_xyz(img.fwd_model, 1);

        G1 = G.Gx;
        G2 = G.Gy;
    otherwise
        error('Invalid dimension');
end

sigma = img.elem_data;
A_matrix = M(img,sigma);

% Jacobi preconditioner - matrix free
d = sqrt(diag(A_matrix));        % vector of diagonal entries

Mfun = @(x) x ./ d;              % left preconditioner  M^{-1} x
Nfun = @(x) x ./ d;              % right preconditioner

% Multiply elementwise by sigma (diagonal) without creating sparse diagonal matrix
% rhs = mu_factor*( dR1 * Sigma * G1 -  dR2* Sigma * G2);
rhs1 = mu_factor*( (dR1{1} .* sigma.')*G1 - (dR2{1} .* sigma.')*G2 );
rhs2 = mu_factor*( (dR1{2} .* sigma.')*G1 - (dR2{2} .* sigma.')*G2 );
rhs3 = mu_factor*( (dR1{3} .* sigma.')*G1 - (dR2{3} .* sigma.')*G2 );

parfor m = 1:num_sensors

    [dlambda_dx(:,m),~,~] = pcg(A_matrix,rhs1(m,:)',1e-9,numel(sigma),Mfun,Nfun);
    [dlambda_dy(:,m),~,~] = pcg(A_matrix,rhs2(m,:)',1e-9,numel(sigma),Mfun,Nfun);
    [dlambda_dz(:,m),~,~] = pcg(A_matrix,rhs3(m,:)',1e-9,numel(sigma),Mfun,Nfun);
end

dlambda = {dlambda_dx,dlambda_dy,dlambda_dz};

end

%% Compute dR
function dR = dRmkj_xyz_optimized(fmdl,j)

% Get model dimension
dim = size(fmdl.nodes,2);

switch dim
    case 2
        error('Implemented only for 3d still!\n');
    case 3
        dR = dRmkj_xyz_optimized_3d(fmdl,j);
    otherwise
        error('Unexpected')
end
end

function dR = dRmkj_xyz(fmdl,j)

% Get model dimension
dim = size(fmdl.nodes,2);

switch dim
    case 2
        error('Implemented only for 3d still!\n');
    case 3
        dR = dRmkj_xyz_3d(fmdl,j);
    otherwise
        error('Unexpected')
end

end

function dR = dRmkj_xyz_optimized_3d(fmdl,j)

persistent ref_pts weights nq

if isempty(ref_pts)
coord = [    0.0000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  1.0000000000000000
    0.0000000000000000  1.0000000000000000  0.0000000000000000
    1.0000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  0.7500000000000000
    0.0000000000000000  0.0000000000000000  0.2500000000000000
    0.0000000000000000  0.7500000000000000  0.0000000000000000
    0.0000000000000000  0.7500000000000000  0.2500000000000000
    0.0000000000000000  0.2500000000000000  0.0000000000000000
    0.0000000000000000  0.2500000000000000  0.7500000000000000
    0.7500000000000000  0.0000000000000000  0.0000000000000000
    0.7500000000000000  0.0000000000000000  0.2500000000000000
    0.7500000000000000  0.2500000000000000  0.0000000000000000
    0.2500000000000000  0.0000000000000000  0.0000000000000000
    0.2500000000000000  0.0000000000000000  0.7500000000000000
    0.2500000000000000  0.7500000000000000  0.0000000000000000
    0.5000000000000000  0.5000000000000000  0.0000000000000000
    0.5000000000000000  0.0000000000000000  0.5000000000000000
    0.5000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.5000000000000000  0.5000000000000000
    0.0000000000000000  0.5000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  0.5000000000000000
    0.2500000000000000  0.2500000000000000  0.0000000000000000
    0.2500000000000000  0.2500000000000000  0.5000000000000000
    0.2500000000000000  0.0000000000000000  0.2500000000000000
    0.2500000000000000  0.0000000000000000  0.5000000000000000
    0.2500000000000000  0.5000000000000000  0.2500000000000000
    0.2500000000000000  0.5000000000000000  0.0000000000000000
    0.0000000000000000  0.2500000000000000  0.2500000000000000
    0.0000000000000000  0.2500000000000000  0.5000000000000000
    0.0000000000000000  0.5000000000000000  0.2500000000000000
    0.5000000000000000  0.2500000000000000  0.2500000000000000
    0.5000000000000000  0.2500000000000000  0.0000000000000000
    0.5000000000000000  0.0000000000000000  0.2500000000000000
    0.2500000000000000  0.2500000000000000  0.2500000000000000];

weights = [  -0.0119047619047619
    -0.0119047619047619
    -0.0119047619047619
    -0.0119047619047619
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.3047619047619048];

ref_pts = coord.';      % 3 × nq
nq = size(ref_pts,2);
end

numElements = size(fmdl.elems,1);
numSensors  = numel(fmdl.sensors);

nodes = fmdl.nodes;
elems = fmdl.elems;

% ---------------------------------------------------------
% Precompute element geometry (ONCE)
% ---------------------------------------------------------
v1  = zeros(3,numElements);
J   = zeros(3,3,numElements);

for k = 1:numElements
    verts = nodes(elems(k,:),:);

    v1(:,k) = verts(1,:).';

    J(:,:,k) = [ ...
        (verts(2,:) - verts(1,:))' , ...
        (verts(3,:) - verts(1,:))' , ...
        (verts(4,:) - verts(1,:))' ];
end

% Vectorized determinant (much faster than loop)
detJ = abs( ...
      J(1,1,:).*(J(2,2,:).*J(3,3,:) - J(2,3,:).*J(3,2,:)) ...
    - J(1,2,:).*(J(2,1,:).*J(3,3,:) - J(2,3,:).*J(3,1,:)) ...
    + J(1,3,:).*(J(2,1,:).*J(3,2,:) - J(2,2,:).*J(3,1,:)) );

detJ = reshape(detJ,numElements,1);

% ---------------------------------------------------------
% Precompute physical quadrature points for ALL elements
% Xi: 3 × nq × numElements
% ---------------------------------------------------------

Xi = zeros(3,nq,numElements);

for k = 1:numElements
    Xi(:,:,k) = v1(:,k) + J(:,:,k)*ref_pts;
end


% ---------------------------------------------------------
% Allocate output
% ---------------------------------------------------------

dRdx = zeros(numSensors,numElements);
dRdy = zeros(numSensors,numElements);
dRdz = zeros(numSensors,numElements);

w = reshape(weights,1,nq,1);   % 1 × nq × 1

sensors = fmdl.sensors;

for m = 1:numSensors

    rm = sensors(m).position';  % 1 x 3

    diff = rm - Xi;              % 3 × nq × numElements

    r2 = sum(diff.^2,1);         % 1 × nq × numElements
    r5 = r2.^(5/2);

    dm_j = diff(j,:,:);

    for p = 1:3

        dm_p = diff(p,:,:);

        num = (j==p)*r2 - 3*dm_j.*dm_p;

        integrand = num ./ r5;

        val = squeeze(sum(integrand .* w,2));  % 1 × numElements

        switch p
            case 1
                dRdx(m,:) = val .* (detJ/6);
            case 2
                dRdy(m,:) = val .* (detJ/6);
            case 3
                dRdz(m,:) = val .* (detJ/6);
        end
    end
end

dR = {dRdx,dRdy,dRdz};

end

function dR = dRmkj_xyz_optimized_2d(fmdl,j)

persistent ref_pts weights nq

if isempty(ref_pts)
coord = [    0.0000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  1.0000000000000000
    0.0000000000000000  1.0000000000000000  0.0000000000000000
    1.0000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  0.7500000000000000
    0.0000000000000000  0.0000000000000000  0.2500000000000000
    0.0000000000000000  0.7500000000000000  0.0000000000000000
    0.0000000000000000  0.7500000000000000  0.2500000000000000
    0.0000000000000000  0.2500000000000000  0.0000000000000000
    0.0000000000000000  0.2500000000000000  0.7500000000000000
    0.7500000000000000  0.0000000000000000  0.0000000000000000
    0.7500000000000000  0.0000000000000000  0.2500000000000000
    0.7500000000000000  0.2500000000000000  0.0000000000000000
    0.2500000000000000  0.0000000000000000  0.0000000000000000
    0.2500000000000000  0.0000000000000000  0.7500000000000000
    0.2500000000000000  0.7500000000000000  0.0000000000000000
    0.5000000000000000  0.5000000000000000  0.0000000000000000
    0.5000000000000000  0.0000000000000000  0.5000000000000000
    0.5000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.5000000000000000  0.5000000000000000
    0.0000000000000000  0.5000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  0.5000000000000000
    0.2500000000000000  0.2500000000000000  0.0000000000000000
    0.2500000000000000  0.2500000000000000  0.5000000000000000
    0.2500000000000000  0.0000000000000000  0.2500000000000000
    0.2500000000000000  0.0000000000000000  0.5000000000000000
    0.2500000000000000  0.5000000000000000  0.2500000000000000
    0.2500000000000000  0.5000000000000000  0.0000000000000000
    0.0000000000000000  0.2500000000000000  0.2500000000000000
    0.0000000000000000  0.2500000000000000  0.5000000000000000
    0.0000000000000000  0.5000000000000000  0.2500000000000000
    0.5000000000000000  0.2500000000000000  0.2500000000000000
    0.5000000000000000  0.2500000000000000  0.0000000000000000
    0.5000000000000000  0.0000000000000000  0.2500000000000000
    0.2500000000000000  0.2500000000000000  0.2500000000000000];

weights = [  -0.0119047619047619
    -0.0119047619047619
    -0.0119047619047619
    -0.0119047619047619
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.3047619047619048];

ref_pts = coord.';      % 3 × nq
nq = size(ref_pts,2);
end

numElements = size(fmdl.elems,1);
numSensors  = numel(fmdl.sensors);

nodes = fmdl.nodes;
elems = fmdl.elems;

% ---------------------------------------------------------
% Precompute element geometry (ONCE)
% ---------------------------------------------------------
v1  = zeros(3,numElements);
J   = zeros(3,3,numElements);

for k = 1:numElements
    verts = nodes(elems(k,:),:);

    v1(:,k) = verts(1,:).';

    J(:,:,k) = [ ...
        (verts(2,:) - verts(1,:))' , ...
        (verts(3,:) - verts(1,:))' , ...
        (verts(4,:) - verts(1,:))' ];
end

% Vectorized determinant (much faster than loop)
detJ = abs( ...
      J(1,1,:).*(J(2,2,:).*J(3,3,:) - J(2,3,:).*J(3,2,:)) ...
    - J(1,2,:).*(J(2,1,:).*J(3,3,:) - J(2,3,:).*J(3,1,:)) ...
    + J(1,3,:).*(J(2,1,:).*J(3,2,:) - J(2,2,:).*J(3,1,:)) );

detJ = reshape(detJ,numElements,1);

% ---------------------------------------------------------
% Precompute physical quadrature points for ALL elements
% Xi: 3 × nq × numElements
% ---------------------------------------------------------

Xi = zeros(3,nq,numElements);

for k = 1:numElements
    Xi(:,:,k) = v1(:,k) + J(:,:,k)*ref_pts;
end


% ---------------------------------------------------------
% Allocate output
% ---------------------------------------------------------

dRdx = zeros(numSensors,numElements);
dRdy = zeros(numSensors,numElements);
dRdz = zeros(numSensors,numElements);

w = reshape(weights,1,nq,1);   % 1 × nq × 1

sensors = fmdl.sensors;

for m = 1:numSensors

    rm = sensors(m).position';  % 1 x 3

    diff = rm - Xi;              % 3 × nq × numElements

    r2 = sum(diff.^2,1);         % 1 × nq × numElements
    r5 = r2.^(5/2);

    dm_j = diff(j,:,:);

    for p = 1:3

        dm_p = diff(p,:,:);

        num = (j==p)*r2 - 3*dm_j.*dm_p;

        integrand = num ./ r5;

        val = squeeze(sum(integrand .* w,2));  % 1 × numElements

        switch p
            case 1
                dRdx(m,:) = val .* (detJ/6);
            case 2
                dRdy(m,:) = val .* (detJ/6);
            case 3
                dRdz(m,:) = val .* (detJ/6);
        end
    end
end

dR = {dRdx,dRdy,dRdz};

end

% THIS IS MISSING OPTIMIZATION
function dR = dRmkj_xyz_3d(fmdl,j)
coord = [    0.0000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  1.0000000000000000
    0.0000000000000000  1.0000000000000000  0.0000000000000000
    1.0000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  0.7500000000000000
    0.0000000000000000  0.0000000000000000  0.2500000000000000
    0.0000000000000000  0.7500000000000000  0.0000000000000000
    0.0000000000000000  0.7500000000000000  0.2500000000000000
    0.0000000000000000  0.2500000000000000  0.0000000000000000
    0.0000000000000000  0.2500000000000000  0.7500000000000000
    0.7500000000000000  0.0000000000000000  0.0000000000000000
    0.7500000000000000  0.0000000000000000  0.2500000000000000
    0.7500000000000000  0.2500000000000000  0.0000000000000000
    0.2500000000000000  0.0000000000000000  0.0000000000000000
    0.2500000000000000  0.0000000000000000  0.7500000000000000
    0.2500000000000000  0.7500000000000000  0.0000000000000000
    0.5000000000000000  0.5000000000000000  0.0000000000000000
    0.5000000000000000  0.0000000000000000  0.5000000000000000
    0.5000000000000000  0.0000000000000000  0.0000000000000000
    0.0000000000000000  0.5000000000000000  0.5000000000000000
    0.0000000000000000  0.5000000000000000  0.0000000000000000
    0.0000000000000000  0.0000000000000000  0.5000000000000000
    0.2500000000000000  0.2500000000000000  0.0000000000000000
    0.2500000000000000  0.2500000000000000  0.5000000000000000
    0.2500000000000000  0.0000000000000000  0.2500000000000000
    0.2500000000000000  0.0000000000000000  0.5000000000000000
    0.2500000000000000  0.5000000000000000  0.2500000000000000
    0.2500000000000000  0.5000000000000000  0.0000000000000000
    0.0000000000000000  0.2500000000000000  0.2500000000000000
    0.0000000000000000  0.2500000000000000  0.5000000000000000
    0.0000000000000000  0.5000000000000000  0.2500000000000000
    0.5000000000000000  0.2500000000000000  0.2500000000000000
    0.5000000000000000  0.2500000000000000  0.0000000000000000
    0.5000000000000000  0.0000000000000000  0.2500000000000000
    0.2500000000000000  0.2500000000000000  0.2500000000000000];

weights = [  -0.0119047619047619
    -0.0119047619047619
    -0.0119047619047619
    -0.0119047619047619
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    -0.0285714285714286
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.0380952380952381
    0.3047619047619048];

numElements = size(fmdl.elems,1);
numSensors  = numel(fmdl.sensors);
numQuadPts  = length(weights);

% --- Precompute element vertices ---
% V: 4 x 3 x numElements
% Stack vertices: 4 x numElements x 3
V = reshape(fmdl.nodes(fmdl.elems',:), [4, numElements, 3]);

v1 = squeeze(V(1,:,:));  % numElements x 3
v2 = squeeze(V(2,:,:));
v3 = squeeze(V(3,:,:));
v4 = squeeze(V(4,:,:));

% Initialize J
J = zeros(3,3,numElements);  % 3 x 3 x numElements

% Fill J using transpose of row vectors
J(:,1,:) = permute(v2 - v1, [2 3 1]);  % 3 x 1 x numElements
J(:,2,:) = permute(v3 - v1, [2 3 1]);
J(:,3,:) = permute(v4 - v1, [2 3 1]);

% Determinants of all J
detJ = zeros(1,numElements);
for k = 1:numElements
    detJ(k) = abs(det(J(:,:,k)));
end

% --- Compute dRdp for all sensors, and for all p ---
dRdx = zeros(numSensors,numElements);
dRdy = zeros(numSensors,numElements);
dRdz = zeros(numSensors,numElements);



for m = 1:numSensors
    rm = fmdl.sensors(m).position;  % 1 x 3

    % Evaluate quadrature points in all elements
    % xi: numQuadPts x 3 x numElements
    xi = reshape(v1', [3,1,numElements]) + ...
        J(:,1,:).*reshape(coord(:,1),[1,numQuadPts,1]) + ...
        J(:,2,:).*reshape(coord(:,2),[1,numQuadPts,1]) + ...
        J(:,3,:).*reshape(coord(:,3),[1,numQuadPts,1]); % 3 x numQuadPts x numElements

    xi = permute(xi,[2 1 3]); % numQuadPts x 3 x numElements

    % dm_vec = rm - xi, same size
    dm_vec = rm - xi; % numQuadPts x 3 x numElements

    dm_norm2 = sum(dm_vec.^2,2);           % numQuadPts x 1 x numElements
    dm_j    = dm_vec(:,j,:);               % numQuadPts x 1 x numElements

    for p = 1:3
        % Compute fun values at all quadrature points
        % ((j==p)*norm(dm)^2 - 3*dm_j*dm_p) / norm(dm)^5

        dm_p    = dm_vec(:,p,:);

        funvals = ((j==p)*dm_norm2 - 3*dm_j.*dm_p) ./ (dm_norm2.^(5/2)); % numQuadPts x 1 x numElements

        % Integrate over quadrature points using weights
        switch p
            case 1
                dRdx(m,:) = squeeze(sum(funvals .* reshape(weights, [numQuadPts,1,1]),1))' .* (detJ/6);
            case 2
                dRdy(m,:) = squeeze(sum(funvals .* reshape(weights, [numQuadPts,1,1]),1))' .* (detJ/6);
            case 3
                dRdz(m,:) = squeeze(sum(funvals .* reshape(weights, [numQuadPts,1,1]),1))' .* (detJ/6);
        end
    end
end

dR = {dRdx,dRdy,dRdz};

end

%% FUNCTIONS: compute_r_matrices_local
function [Rx,Ry,Rz,fmdl] = compute_r_matrices_local(fmdl,sensor_locations)

% Get model dimension
dim = size(fmdl.nodes,2);

switch dim
    case 2
        error('Implemented only for 3d still!\n');
    case 3
        [Rx,Ry,Rz,fmdl]  = compute_r_matrices_local_3d(fmdl,sensor_locations);
    otherwise
        error('Unexpected')
end
end

function [Rx,Ry,Rz,fmdl] = compute_r_matrices_local_3d(fmdl,sensor_locations)

numElements = size(fmdl.elems,1);
numSensors = size(sensor_locations,1);

% ------------------------------------------------------------
% Persistent quadrature rule (avoid reallocation every call)
% ------------------------------------------------------------
persistent ref_pts weights nq;

if isempty(ref_pts)

    coord = [    0.0000000000000000  0.0000000000000000  0.0000000000000000
        0.0000000000000000  0.0000000000000000  1.0000000000000000
        0.0000000000000000  1.0000000000000000  0.0000000000000000
        1.0000000000000000  0.0000000000000000  0.0000000000000000
        0.0000000000000000  0.0000000000000000  0.7500000000000000
        0.0000000000000000  0.0000000000000000  0.2500000000000000
        0.0000000000000000  0.7500000000000000  0.0000000000000000
        0.0000000000000000  0.7500000000000000  0.2500000000000000
        0.0000000000000000  0.2500000000000000  0.0000000000000000
        0.0000000000000000  0.2500000000000000  0.7500000000000000
        0.7500000000000000  0.0000000000000000  0.0000000000000000
        0.7500000000000000  0.0000000000000000  0.2500000000000000
        0.7500000000000000  0.2500000000000000  0.0000000000000000
        0.2500000000000000  0.0000000000000000  0.0000000000000000
        0.2500000000000000  0.0000000000000000  0.7500000000000000
        0.2500000000000000  0.7500000000000000  0.0000000000000000
        0.5000000000000000  0.5000000000000000  0.0000000000000000
        0.5000000000000000  0.0000000000000000  0.5000000000000000
        0.5000000000000000  0.0000000000000000  0.0000000000000000
        0.0000000000000000  0.5000000000000000  0.5000000000000000
        0.0000000000000000  0.5000000000000000  0.0000000000000000
        0.0000000000000000  0.0000000000000000  0.5000000000000000
        0.2500000000000000  0.2500000000000000  0.0000000000000000
        0.2500000000000000  0.2500000000000000  0.5000000000000000
        0.2500000000000000  0.0000000000000000  0.2500000000000000
        0.2500000000000000  0.0000000000000000  0.5000000000000000
        0.2500000000000000  0.5000000000000000  0.2500000000000000
        0.2500000000000000  0.5000000000000000  0.0000000000000000
        0.0000000000000000  0.2500000000000000  0.2500000000000000
        0.0000000000000000  0.2500000000000000  0.5000000000000000
        0.0000000000000000  0.5000000000000000  0.2500000000000000
        0.5000000000000000  0.2500000000000000  0.2500000000000000
        0.5000000000000000  0.2500000000000000  0.0000000000000000
        0.5000000000000000  0.0000000000000000  0.2500000000000000
        0.2500000000000000  0.2500000000000000  0.2500000000000000];

    weights = [  -0.0119047619047619
        -0.0119047619047619
        -0.0119047619047619
        -0.0119047619047619
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        -0.0285714285714286
        -0.0285714285714286
        -0.0285714285714286
        -0.0285714285714286
        -0.0285714285714286
        -0.0285714285714286
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.0380952380952381
        0.3047619047619048];
    ref_pts = coord.';      % 3 × nq
    nq = size(ref_pts,2);
end

nodes = fmdl.nodes;
elems = fmdl.elems;

% ------------------------------------------------------------
% Precompute element geometry (done once!)
% ------------------------------------------------------------

v1  = zeros(3,numElements);
Jall = zeros(3,3,numElements);
detJ = zeros(numElements,1);

for k = 1:numElements

    verts = nodes(elems(k,:),:);

    v1(:,k) = verts(1,:).';

    J = [ ...
        (verts(2,:) - verts(1,:))' , ...
        (verts(3,:) - verts(1,:))' , ...
        (verts(4,:) - verts(1,:))' ];

    Jall(:,:,k) = J;
    detJ(k)     = abs(det(J));
end

% ------------------------------------------------------------
% Allocate outputs
% ------------------------------------------------------------

Rx = zeros(numSensors,numElements);
Ry = zeros(numSensors,numElements);
Rz = zeros(numSensors,numElements);

% ------------------------------------------------------------
% Main loop over sensors (parallelizable)
% ------------------------------------------------------------
parfor m = 1:numSensors   % <-- change to "for" if no Parallel Toolbox

    rm = sensor_locations(m,:).';

    Rx_local = zeros(1,numElements);
    Ry_local = zeros(1,numElements);
    Rz_local = zeros(1,numElements);

    for k = 1:numElements

        % Map quadrature points to physical tetrahedron
        Xi = v1(:,k) + Jall(:,:,k) * ref_pts;    % 3 × nq

        diff = rm - Xi;                         % 3 × nq
        r2   = sum(diff.^2,1);                  % 1 × nq

        kernel = diff ./ (r2.^1.5);             % 3 × nq

        R = (detJ(k)/6) * (kernel * weights);  % 3 × 1

        Rx_local(k) = R(1);
        Ry_local(k) = R(2);
        Rz_local(k) = R(3);
    end

    Rx(m,:) = Rx_local;
    Ry(m,:) = Ry_local;
    Rz(m,:) = Rz_local;
end

% ------------------------------------------------------------
% Store
% ------------------------------------------------------------

fmdl.R.Rx = Rx;
fmdl.R.Ry = Ry;
fmdl.R.Rz = Rz;

end


%% FUNCTIONS: compute_gamma_matrices_local
function img = compute_gamma_matrices_local(img)

mu_factor = img.fwd_model.mu0/(4*pi);

num_sensors = numel(img.fwd_model.sensors);

sensor_locations = zeros(numel(img.fwd_model.sensors),3);

for i = 1: numel(img.fwd_model.sensors)
    sensor_locations(i,:) = img.fwd_model.sensors(i).position;
end

[Rx,Ry,Rz,fmdl] = compute_r_matrices_local(img.fwd_model,sensor_locations);

% THIS IS CRUCIAL!!!!! so the R matrices are updated whenever compute_gamma_matrices_local is called, otherwise they will be the same as the initial img.fwd_model, which might be wrong!
img.fwd_model = fmdl;

% Convenience handles
R.Rx = Rx;
R.Ry = Ry;
R.Rz = Rz;
G = img.fwd_model.G;

% Sigma = sparse(1:length(img.elem_data), 1:length(img.elem_data), img.elem_data);
Sigma = spdiags(img.elem_data(:), 0, length(img.elem_data), length(img.elem_data));

% NEW: The matrix g_{dl}^m, gives the components of the measurement axis of
% sensor m on the canonical R^3 basis

g = zeros(num_sensors,3,3);
for m = 1:numel(img.fwd_model.sensors)
    g(m,:,:) = [...
        img.fwd_model.sensors(m).axes.axis1;
        img.fwd_model.sensors(m).axes.axis2;
        img.fwd_model.sensors(m).axes.axis3];
end

Cx = ( -R.Rz * Sigma * G.Gy +  R.Ry * Sigma * G.Gz );
Cy = ( -R.Rx * Sigma * G.Gz +  R.Rz * Sigma * G.Gx );
Cz = ( -R.Ry * Sigma * G.Gx +  R.Rx * Sigma * G.Gy );

Gamma1 = mu_factor*(g(:,1,1).*Cx + g(:,1,2).*Cy + g(:,1,3).*Cz);
Gamma2 = mu_factor*(g(:,2,1).*Cx + g(:,2,2).*Cy + g(:,2,3).*Cz);
Gamma3 = mu_factor*(g(:,3,1).*Cx + g(:,3,2).*Cy + g(:,3,3).*Cz);

img.Gamma1 = Gamma1;
img.Gamma2 = Gamma2;
img.Gamma3 = Gamma3;

end

%% FUNCTIONS: calc_jacobian_mdeit_local
function [J,img] = calc_jacobian_mdeit_local(img,A,recon_mode,select_sensor_axis)

if nargin<4
    %Check if mode is mdeit3
    assert(strcmp(recon_mode,'mdeit3'));
end


valid_recon_modes = {'mdeit1','mdeit3'};
assert(ismember(lower(recon_mode),valid_recon_modes));

switch recon_mode

    case 'mdeit1'
        J = calc_jacobian_1axis_local(img,A,select_sensor_axis);

    case 'mdeit3'
        [Jx,Jy,Jz] = calc_jacobian_3axis_local(img,A);
        J = [Jx;Jy;Jz];
end

end

function J = calc_jacobian_1axis_local(img,A,select_sensor_axis)

mu0 = img.fwd_model.mu0;

n_nodes = size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);

% Compute Gamma matrices
img = compute_gamma_matrices_local(img);

R = img.fwd_model.R;
G = img.fwd_model.G;

switch select_sensor_axis
    case 1
        Gamma = img.Gamma1;
        R1 = R.Rz.';
        R2 = R.Ry.';
    case 2
        Gamma = img.Gamma2;
        R1 = R.Rx.';
        R2 = R.Rz.';
    case 3
        Gamma = img.Gamma3;
        R1 = R.Ry.';
        R2 = R.Rx.';
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

function [Jx,Jy,Jz] = calc_jacobian_3axis_local(img,A)

mu0 = img.fwd_model.mu0;
n_nodes = size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);

% Compute Gamma matrices
img = compute_gamma_matrices_local(img);

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



function sensor_locations = fetch_sensor_locations(img)
n_sensors = numel(img.fwd_model.sensors);
sensor_locations = zeros(n_sensors,3);

for m = 1: numel(img.fwd_model.sensors)
    sensor_locations(m,:) = img.fwd_model.sensors(m).position;
end

end
