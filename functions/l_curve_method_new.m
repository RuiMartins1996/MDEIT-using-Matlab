function [lambda_opt,dx] = l_curve_method_new(A,r,lambda_vector,U,S,V)

% Compute SVD if not provided
if nargin<4
    fprintf('Computing SVD\n');
    [U,S,V] = svd(A,'econ');
    fprintf('Done\n');
end

%Filter factors (several things can be expressed as a function of them)
f = @(lambda) 1./(1+lambda./(diag(S).^2));
df = @(lambda) 1/lambda .*(f(lambda)-1).*f(lambda);

% Compute residual norms and x norms for all lambda
rho = zeros(1,length(lambda_vector));
eta = zeros(1,length(lambda_vector));

beta = U.'*r;
sigma = diag(S);

for i = 1:length(lambda_vector)

    lambda = lambda_vector(i);

    % x_norms(i) = norm(d_sigma,2);
    rho(i) = sum(((1-f(lambda)).*beta).^2) ;% eq(9) Hansen
    eta(i) = sum((f(lambda).*beta./sigma).^2); % eq(8) Hansen
end

%% Compute optimal hyperparameter 

% Sort by residual norm (important for monotone x)
[rho_sorted, idx] = sort(rho);
eta_sorted = eta(idx);

% Log-log scale for L-curve
lrho = log(rho_sorted);
leta = log(eta_sorted);

lambda_vector_dense = logspace(log10(min(lambda_vector)), log10(max(lambda_vector)), 100);

leta_dense = interp1(lambda_vector, leta, lambda_vector_dense, 'pchip');
lrho_dense = interp1(lambda_vector, lrho, lambda_vector_dense, 'pchip');
leta_dense = leta_dense(:);
lrho_dense = lrho_dense(:);

rho_dense = exp(lrho_dense);
eta_dense = exp(leta_dense);

F = @(lambda) spdiag(f(lambda));
deta_dense_func = @(lambda) (2/lambda) .* ...
    (U.'*r).'*...
    (F(lambda)-speye(numel(diag(S))))*...
    F(lambda).^2*....
    spdiag(1./diag(S)).^2*...
    U.'*r;

deta_dense = zeros(numel(lambda_vector_dense),1);
for n = 1:length(lambda_vector_dense)
    deta_dense(n) = deta_dense_func(lambda_vector_dense(n));
end

kappa_dense = -2.*(rho_dense.*eta_dense).*...
    (...
    lambda_vector_dense(:).^2.*eta_dense+...
    lambda_vector_dense(:).*rho_dense+...
    rho_dense.*eta_dense./deta_dense)./...
    ((lambda_vector_dense(:).^2.*eta_dense.^2+rho_dense.^2).^(3/2));

% Find maximum curvature (corner of L-curve)
[kappa_max, imax] = max(kappa_dense);
lambda_opt = lambda_vector_dense(imax);
opt_r = rho_dense(imax);
opt_x = eta_dense(imax);

fprintf('Optimal residual norm = %.4e\n', opt_r);
fprintf('Optimal solution norm = %.4e\n', opt_x);
fprintf('Maximum curvature = %.4e\n', kappa_dense(imax));
fprintf('Optimal hyperparameter = %.4e\n', lambda_opt);


%% Reconstruct

dx = V * ((sigma ./ (sigma.^2 + lambda_opt)) .* U.' * r);


%% Plot
figure
subplot(1,2,1)
hold on
plot(rho_sorted,eta_sorted,'b.');
plot(rho,eta);
plot(opt_r,opt_x,'r.','MarkerSize',10)
text(opt_r+0.1,opt_x,strcat('$\lambda = ',num2str(lambda_opt),'$'),'Interpreter','latex')
hold off
grid on;grid minor;box on;
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('$\rho$','Interpreter','latex')
ylabel('$\eta$','Interpreter','latex')

subplot(1,2,2)
plot(lambda_vector_dense,kappa_dense,'.')
set(gca,'XScale','log')
grid on;grid minor;box on;
xlabel('$\lambda$','Interpreter','latex')
ylabel('$\kappa(\lambda)$','Interpreter','latex')


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