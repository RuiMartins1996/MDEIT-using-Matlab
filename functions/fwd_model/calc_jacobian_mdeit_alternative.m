function [J,img] = calc_jacobian_mdeit_alternative(img,recon_mode,select_sensor_axis)

check_input(recon_mode);

if isfield(img,'jacobian')
    J = img.jacobian;
    return;
end

switch recon_mode
    case 'eit'
        J = calc_jacobian_eit(img);
    case 'mdeit1'
        J = calc_jacobian_1axis(img,select_sensor_axis);
    case 'mdeit3'
        error('Unimplemented');
end
end


function check_input(recon_mode)

% Check if recon_mode is valid
valid_modes = {'mdeit1', 'mdeit3','eit'};
if ~ismember(recon_mode, valid_modes)
    error('invalid recon_mode "%s". Must be''mdeit1'', or ''mdeit3''.', recon_mode);
end

end


% These versions are optimized with respect to former ones (avoids
% permutation and solves the adjoint systems with \, not pcg)
function J = calc_jacobian_1axis(img,select_sensor_axis)

mu0 = img.fwd_model.mu0;

n_nodes = size(img.fwd_model.nodes,1);
n_elem = size(img.fwd_model.elems,1);

num_stim = numel(img.fwd_model.stimulation);
num_sensors = numel(img.fwd_model.sensors);
num_electrodes = numel(img.fwd_model.electrode);

J = zeros(num_stim*num_sensors,n_elem);

[I1,I2,I3] = compute_w_matrices(img.fwd_model,mu0,select_sensor_axis);
[G1,G2,G3,element_volume] = compute_gradient_matrix(img.fwd_model);

img = compute_gamma_matrices(img);

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

P = lambda.';

for s = 1:num_stim
    Js = ...
        (I1+P*(G1.')*diag(element_volume))*diag(G1*u(:,s)) + ...
        (I2+P*(G2.')*diag(element_volume))*diag(G2*u(:,s)) + ...
        (I3+P*(G3.')*diag(element_volume))*diag(G3*u(:,s));

    J(num_sensors*(s-1)+1:s*num_sensors,:) = Js;
end 

return
end






%% Functions
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