function [mu_min,optimal_id,V_lambda,dx] = generalized_cross_validation(imdl,data,lambda_vector,opts)
%generalized_cross_validation: Obtain the optimal regularization parameter
%through the generalized cross validation (GCV) method described in "Generalized Cross-Validation
% as a Method for Choosing a Good Ridge Parameter"
%
%   Inputs:
%       imdl - MDEIT inverse model
%       
%       data - Difference data for MDEIT (difference of magnetic fields between
%       inhomogeneous and homogeneous forward models)
%       
%       mu_vector - Vector of regularization parameters at which to compute
%       the cost function of GCV. Notice that the forward model assumed in
%       the GCV article considers n\lambda as the regularization, where n
%       is the number of rows of the Jacobian matrix. mu_vector gives an
%       array of values of \lambda.
%       
%       options - struct() that controls the behaviour of this function:
%           -> options.reconstruct (boolean): Controls wether the we should
%           perform conductivity reconstruction with the optimal
%           regularization parameter in this function (default = false)
%
%
%   Outputs:
%       mu_min - Optimal regularization parameter
%       
%       optimal_id - Index of mu_vector/V_mu at which the minimum is located       
%       
%       V_mu - minimum of the GCV cost function
%       
%       dx - Conductivity difference between inhomogeneous and homogeneous model

arguments
    imdl struct
    data (:,1) double % column vector of doubles
    lambda_vector (:,1) double % column vector of doubles
    opts.reconstruct logical = false;
end

if isfield(imdl,'recon_type')
    assert(strcmp(imdl.recon_type,'difference'),'imdl.recon_type must be set to "difference"');
elseif isfield(imdl,'reconst_type')
    assert(strcmp(imdl.reconst_type,'difference'),'imdl.reconst_type must be set to "difference"');
else
    error('imdl.reconst_type or imdl.recon_type must be set to "difference"');
end
    

assert(isfield(imdl,'jacobian_bkgnd'),'Must have the field imdl.jacobian_bkgnd');

jacobian_bkgnd = imdl.jacobian_bkgnd;
n_elems = size(imdl.fwd_model.elems,1);

assert(isfield(imdl,'recon_mode'),'Expected field imdl.recon_mode');
recon_mode = imdl.recon_mode;

valid_modes = {'mdeit1', 'mdeit3','eit'};
assert(ismember(recon_mode, valid_modes),'my_inv_solve: invalid recon_mode "%s". Must be ''eit'',''mdeit1'' or ''mdeit3''.', recon_mode)

n_sensors = numel(imdl.fwd_model.sensors);
n_inj = numel(imdl.fwd_model.stimulation);


if strcmp(recon_mode,'mdeit1')
    assert(isfield(imdl,'select_sensor_axis'),'Expected field imdl.select_sensor_axis');
    select_sensor_axis = imdl.select_sensor_axis;

    % Check if data size matches expected size:
    assert(size(data,1) == n_sensors*n_inj,'Size of data is wrong');
end

if strcmp(recon_mode,'mdeit3')
    % Check if data size matches expected size:
    assert(size(data,1) == 3*n_sensors*n_inj,'Size of data is wrong');
end

if strcmp(recon_mode,'eit')
    % Check if data size matches expected size:
    n_meas = 0;

    for i = 1:numel(imdl.fwd_model.stimulation)
        n_meas = n_meas +  size(imdl.fwd_model.stimulation(i).meas_pattern,1);
    end
    assert(size(data,1) == n_meas,'Size of data is wrong');
end

% Make background image
if strcmp(recon_mode,'mdeit1') || strcmp(recon_mode,'mdeit3')
    imgh = mk_image_mdeit(imdl.fwd_model,jacobian_bkgnd.value);
elseif strcmp(recon_mode,'eit')
    imgh = mk_image(imdl.fwd_model,jacobian_bkgnd.value);
else
    error('Should not have reached this branch');
end

if isstruct(jacobian_bkgnd) && isfield(jacobian_bkgnd,'value')
    x0 = jacobian_bkgnd.value*ones(n_elems,1);
elseif isnumeric(jacobian_bkgnd) && isvector(jacobian_bkgnd)
    x0 = jacobian_bkgnd;
else
    error('imdl.jacobian_bkgnd has unexpectd format');
end

fprintf('Running GCV: %s\n',recon_mode);

% Compute linearized jacobian
fprintf('Computing Jacobian\n');

lambdatimesdAdp = @(lambda) computeLambdaTimesDaDp(imgh,lambda);
A = @(sigma) M(imgh,sigma);

if strcmp(recon_mode,'mdeit1')
    jac = @(x) calc_jacobian_mdeit(...
        imgh, x,lambdatimesdAdp,A,recon_mode,select_sensor_axis);
    J = jac(x0);
elseif strcmp(recon_mode,'mdeit3')
    jac = @(x) calc_jacobian_mdeit(imgh, x,lambdatimesdAdp,A,recon_mode);
    J = jac(x0);
elseif strcmp(recon_mode,'eit')
    J = calc_jacobian(imgh);
end

% The GCV function can be computed much faster in terms of the SVD of J
% without inverting J'*J+n*lambda*speye(n_elems) for every single lambda

% Check Generalized Cross-Validation as a Method for Choosing a Good Ridge
% Parameter for details

fprintf('Precomputing SVD...\n')

tic 
[U,S,V] = svd(J,'econ');

fprintf('Time %2.2f\n',toc);
sigma = diag(S);             % singular values
Uy = U' * data;              % coordinates of data in U-basis
m = size(J,1);
n = length(lambda_vector);

V_lambda = zeros(n,1);

for i = 1:n
    lambda = lambda_vector(i);

    gamma = sigma.^2 ./ (sigma.^2 + m*lambda);    % shrinkage factors

    one_minus_gamma = 1 - gamma;

    numerator = (1/m) * sum( (one_minus_gamma .* Uy).^2 );

    denominator = ( (1/m) * sum(one_minus_gamma) )^2;

    V_lambda(i) = numerator / denominator;
end

optimal_id = find(V_lambda == min(V_lambda));

% V1 = V_lambda;

% TRY randomized-SVD algorithm (THIS APPROACH ENDS UP TAKING LONGER THAN
% JUST svds(J,'econ');
% tic
% Compute the first singular value (rough estimate, to avoid using svds(J,1))
% niter = 10;
% sigma1 = estimate_sigma1(J, niter);
% sigma1 = normest(J);

% Set tolerance to 1% of maximum singular value
% [U,S,V] = rsvd(J,[],0.01*sigma1);

% [U,S,V] = rsvd(J,m);
% 
% fprintf('Time %2.2f\n',toc);
% 
% sigma = diag(S);             % singular values
% Uy = U' * data;              % coordinates of data in U-basis
% m = size(J,1);
% n = length(lambda_vector);
% 
% V_lambda = zeros(n,1);
% 
% for i = 1:n
%     lambda = lambda_vector(i);
% 
%     gamma = sigma.^2 ./ (sigma.^2 + m*lambda);    % shrinkage factors
% 
%     one_minus_gamma = 1 - gamma;
% 
%     numerator = (1/m) * sum( (one_minus_gamma .* Uy).^2 );
% 
%     denominator = ( (1/m) * sum(one_minus_gamma) )^2;
% 
%     V_lambda(i) = numerator / denominator;
% end
% 
% V2 = V_lambda;
% 
% rank(J)
% optimal_id = find(V_lambda == min(V_lambda));

% TRY gpuArray (THIS APPROACH ENDS UP TAKING LONGER TOO)
% tic
% Jg = gpuArray(J);
% [U,S,V] = svd(Jg,'econ');
% fprintf('Time %2.2f\n',toc);

if numel(optimal_id)>1
    warning('Multiple minimum values found! Taking smallest one');
    optimal_id = optimal_id(V_lambda(optimal_id) == min(V_lambda(optimal_id)));
end

% cla
% hold on
% plot(V_lambda)
% plot(V_lambda(optimal_id),'r.','MarkerSize',10)
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% hold off

id_min = find(lambda_vector == min(lambda_vector));
id_max = find(lambda_vector == max(lambda_vector));

if optimal_id == id_min || optimal_id == id_max
    warning('Not a local minimum!');
end

mu_min = m*lambda_vector(optimal_id);

if opts.reconstruct == true
    % lhs = J'*J+mu_min*speye(n_elems); rhs = -J'*data;
    
    fprintf('Computing inverse solution for optimal reg. parameter \n');
    dx = -V * ((sigma ./ (sigma.^2 + mu_min)) .* Uy);

else
    dx = []; % Explicitly set to empty
end

end

%% FUNCTIONS
function out = M(img,sigma)

numNodes = size(img.fwd_model.nodes,1);

img.elem_data = sigma;
s_mat = system_mat_1st_order(img);

Ac = s_mat.E(1:numNodes,1:numNodes);
Ae = s_mat.E(1:numNodes,numNodes+1:end);
Ad = s_mat.E(numNodes+1:end,numNodes+1:end);

out = Ac-Ae*inv(Ad)*Ae';
end


function sigma1 = estimate_sigma1(A, niter)
    if nargin < 2, niter = 20; end
    x = randn(size(A,2),1);
    for i = 1:niter
        x = A'*(A*x);        % multiply by A'A
        x = x / norm(x);
    end
    sigma1 = norm(A*x);      % ||A*x|| ≈ largest singular value
end