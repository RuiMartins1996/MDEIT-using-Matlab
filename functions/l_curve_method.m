function s = l_curve_method(imgh,imgi,recon_mode,lambda_vec)

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
    lambda_vec = logspace(-2,1,5);
else
    assert(all(lambda_vec>0),'lambda_vec must be positive')
end

fprintf('Running L-curve method ...\n')

% Allocate for residual and solution norms
residual_norms = zeros(1,length(lambda_vec));
x_norms = zeros(1,length(lambda_vec));

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

% Define residual ||J(sigma_0)*Delta_sigma - y||_2
lambdatimesdAdp = @(x) computeLambdaTimesDaDp(imgi,x);
A = @(x) M(imgi,x);

jac = @(x) calc_jacobian_mdeit(imgi,x,lambdatimesdAdp,A,recon_mode,select_sensor_axis);
J = jac(imgh.elem_data);

res = @(x,y) norm(J*x-y,2); 

% Make inverse model
imdl= eidors_obj('inv_model','my_inverse_model');

imdl.fwd_model = imgi.fwd_model;
imdl.jacobian_bkgnd = struct('value',1.0);
imdl.solve = @eidors_default;
imdl.RtR_prior = @(x,Jx) speye(numel(x));
imdl.max_iterations = 1;
imdl.recon_mode = recon_mode;


for i = 1:length(lambda_vec)
    
    imdl.hyperparameter = struct('value',lambda_vec(i));

    img_output = inv_solve_mdeit(imdl,datah,datai,J);
    data_diff = datai-datah;
    d_sigma = img_output.elem_data;
    
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

% --- Option: use log-log scale if doing L-curve analysis ---
lr = log(rs);
lx = log(xs);

% Fit a cubic smoothing spline: p controls smoothness (0–1)
p = 0.9;  % try values 0.8–0.99 depending on noise
pp = csaps(lr, lx, p);

% Evaluate on a dense grid
lr_dense = linspace(min(lr), max(lr), 200);
lx_dense = fnval(pp, lr_dense);

% Derivatives of the spline
d1 = fnval(fnder(pp,1), lr_dense);
d2 = fnval(fnder(pp,2), lr_dense);

% Curvature formula: κ = |y''| / (1 + y'^2)^(3/2)
kappa = abs(d2) ./ (1 + d1.^2).^(3/2);

% Interpolate curvature values from lr_dense to the original lr points
kappa_interp = interp1(lr_dense, kappa, lr, 'pchip');

% Find maximum curvature (corner of L-curve)
[~, imax] = max(kappa_interp);
opt_lr = lr_dense(imax);
opt_lx = lx_dense(imax);
lambda_opt = lambda_vec(imax);

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

% function ltdAdsigma = computeLambdaTimesDaDp(img,lambda)
% 
% numNodes = size(img.fwd_model.nodes,1);
% numElements = size(img.fwd_model.elems,1);
% 
% numNodesPerElement = size(img.fwd_model.elems(1,:),2);
% idj = zeros(numNodesPerElement*numElements,1);
% idk = zeros(numNodesPerElement*numElements,1);
% values = zeros(numNodesPerElement*numElements,1);
% 
% for k = 1:numElements
%     nodeIds = img.fwd_model.elems(k,:);
% 
%     %No need to recompute the gradients at this element, they are stored in
%     %the G matrices
% 
%     Gk = diag(lambda(nodeIds))*(...
%         img.fwd_model.G.Gx(k,nodeIds)'*img.fwd_model.G.Gx(k,nodeIds)+...
%         img.fwd_model.G.Gy(k,nodeIds)'*img.fwd_model.G.Gy(k,nodeIds)+...
%         img.fwd_model.G.Gz(k,nodeIds)'*img.fwd_model.G.Gz(k,nodeIds));
% 
%     ids = numNodesPerElement*(k-1)+1:numNodesPerElement*k;
%     idj(ids) = nodeIds;
%     idk(ids) = k;
%     % dont forget to multiply by volume
%     values(ids) = img.fwd_model.elem_volume(k)*sum(Gk,1);
% end
% 
% ltdAdsigma = sparse(idj,idk,values,numNodes,numElements);
% end

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