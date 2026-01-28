function [data,u] = fwd_solve_mdeit(img,u)

if img.fwd_solve.get_all_meas ~=1
    error('get_all_meas should be set to 1 for u.volt field to exist')
end

if nargin<2
    %% Forward solve the EIT model
    u = fwd_solve(img);
end

%% Output
img = compute_gamma_matrices(img);

Bx = img.Gamma1*u.volt;
By = img.Gamma2*u.volt;
Bz = img.Gamma3*u.volt;

data = struct('Bx',Bx,'By',By,'Bz',Bz);
end

