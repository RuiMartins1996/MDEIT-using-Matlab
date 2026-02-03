function s = l_curve_method(fmdl_reconstruction,imgh,imgi,recon_mode,lambda_vector)

valid_modes = {'eit', 'mdeit1', 'mdeit3'};
if ~ismember(recon_mode, valid_modes)
    error('invalid recon_mode "%s". Must be ''eit'', ''mdeit1'', or ''mdeit3''.', recon_mode);
end

if strcmp(recon_mode,'mdeit1')
    select_sensor_axis = 1;
else
    select_sensor_axis = [];
end

if nargin<3
    lambda_vector = logspace(-2,1,5);
else
    assert(all(lambda_vector>0),'lambda_vec must be positive')
end

fprintf('Running L-curve method ...\n')

% Allocate for residual and solution norms
residual_norms = zeros(1,length(lambda_vector));
x_norms = zeros(1,length(lambda_vector));

% Forward solve img to get data
if strcmp(recon_mode,'mdeit1')
    datah = fwd_solve_mdeit(imgh);
    datai = fwd_solve_mdeit(imgi);

    datah = datah.Bx(:);
    datai = datai.Bx(:);
elseif strcmp(recon_mode,'mdeit3')
    
    datah = fwd_solve_mdeit(imgh);
    datai = fwd_solve_mdeit(imgi);

    datah = [datah.Bx(:);datah.By(:);datah.Bz(:)];
    datai = [datai.Bx(:);datai.By(:);datai.Bz(:)];    
elseif strcmp(recon_mode,'eit')
    error('Not implemented yet')
end

img_reconstruction = mk_image_mdeit(fmdl_reconstruction,(max(imgh.elem_data)));

% Define residual ||J(sigma_0)*Delta_sigma - y||_2
lambdatimesdAdp = @(x) computeLambdaTimesDaDp(img_reconstruction,x);
A = @(x) M(img_reconstruction,x);

jac = @(x) calc_jacobian_mdeit(img_reconstruction,x,lambdatimesdAdp,A,recon_mode,select_sensor_axis);
J = jac(img_reconstruction.elem_data);

res = @(x,y) norm(J*x-y,2); 

% Make inverse model
imdl= eidors_obj('inv_model','my_inverse_model');

imdl.fwd_model = fmdl_reconstruction;
imdl.jacobian_bkgnd = struct('value',1.0);
imdl.solve = @eidors_default;
imdl.RtR_prior = @(x,Jx) x;
imdl.recon_mode = recon_mode;
imdl.recon_type = "difference";
imdl.jacobian_bkgnd = struct('value',max(imgh.elem_data));

for i = 1:length(lambda_vector)
    
    imdl.hyperparameter = struct('value',lambda_vector(i));

    img_output = inv_solve_mdeit(imdl,datah,datai,J);
    data_diff = datai-datah;
    d_sigma = img_output.elem_data;
    
    figure('Name',strcat(recon_mode,'| lambda = ',num2str(lambda_vector(i))));
    show_fem(img_output);
    pause(1e-10)

    residual_norms(i) = res(d_sigma,data_diff);
    x_norms(i) = norm(img_output.elem_data,2);
end

%% TODO : Compute the curvature according to L-curve article

%% Compute the curvature by fitting a smoothing cubic spline to graph

% Given data
r = residual_norms(:);
xnorm = x_norms(:);

% Sort by residual norm (important for monotone x)
[rs, idx] = sort(r);
xs = xnorm(idx);

% Log-log scale for L-curve
lr = log(rs);
lx = log(xs);

% --- Use monotone-preserving piecewise cubic interpolation ---
lr_dense = linspace(min(lr), max(lr), 200);
lx_dense = interp1(lr, lx, lr_dense, 'pchip');

% --- Numerical derivatives using finite differences ---
d1 = gradient(lx_dense, lr_dense);          % first derivative
d2 = gradient(d1, lr_dense);               % second derivative

% --- Curvature formula: κ = |y''| / (1 + y'^2)^(3/2)
kappa = abs(d2) ./ (1 + d1.^2).^(3/2);

% --- Interpolate curvature back to original lr points if needed
kappa_interp = interp1(lr_dense, kappa, lr, 'pchip');

% Find maximum curvature (corner of L-curve)
[~, imax] = max(kappa_interp);
opt_lr = lr_dense(imax);
opt_lx = lx_dense(imax);

% Assuming lambda_vector corresponds to lr_dense
lambda_opt = lambda_vector(imax);

% Convert back to linear scale
opt_r = exp(opt_lr);
opt_x = exp(opt_lx);

fprintf('Optimal residual norm = %.4e\n', opt_r);
fprintf('Optimal solution norm = %.4e\n', opt_x);
fprintf('Maximum curvature = %.4e\n', kappa_interp(imax));
fprintf('Optimal hyperparameter = %.4e\n', lambda_opt);


% Output

s = struct();

s.lambda_opt = lambda_opt;
s.residual_norms = residual_norms;
s.x_norms = x_norms;
s.optimal_id = imax;

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