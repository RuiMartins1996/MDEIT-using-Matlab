function [m,u] = fwd_solve_mdeit_alternative(img,u)

if img.fwd_solve.get_all_meas ~=1
    error('get_all_meas should be set to 1 for u.volt field to exist')
end

if nargin<2
    %% Forward solve the EIT model
    u = fwd_solve(img);
end

n_elem = numel(img.elem_data);
n_sensors = numel(img.fwd_model.sensors);
n_stim = numel(img.fwd_model.stimulation);

[I1,I2,I3] = compute_w_matrices(img.fwd_model,img.fwd_model.mu0,3);

[G1,G2,G3] = compute_gradient_matrix(img.fwd_model);

Sigma = sparse(1:n_elem,1:n_elem,img.elem_data,n_elem,n_elem);

m = (I1*Sigma*G1 + I2*Sigma*G2 + I3*Sigma*G3)*u.volt;

end

