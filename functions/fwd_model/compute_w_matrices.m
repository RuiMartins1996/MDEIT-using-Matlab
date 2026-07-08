function [I1,I2,I3] = compute_w_matrices(fmdl,mu0,select_sensor_axis)

assert(isfield(fmdl,'sensors'));
assert(ismember(select_sensor_axis,[1,2,3]))

n_sensors = numel(fmdl.sensors);

% Compute R matrices with default quadrature rule
[Rx,Ry,Rz] = compute_r_matrices(fmdl);

% Fetch the measurement axis from the sensors

nx_diag = zeros(n_sensors,1);
ny_diag = zeros(n_sensors,1);
nz_diag = zeros(n_sensors,1);

for m = 1:n_sensors
    sensor = fmdl.sensors(m);

    switch select_sensor_axis
        case 1
            sensor_axis = sensor.axes.axis1;
        case 2
            sensor_axis = sensor.axes.axis2;
        case 3
            sensor_axis = sensor.axes.axis3;
    end

    nx_diag(m) = dot(sensor_axis,[1,0,0]);
    ny_diag(m) = dot(sensor_axis,[0,1,0]);
    nz_diag(m) = dot(sensor_axis,[0,0,1]);
end

Nx = sparse(1:n_sensors,1:n_sensors,nx_diag,n_sensors,n_sensors);
Ny = sparse(1:n_sensors,1:n_sensors,ny_diag,n_sensors,n_sensors);
Nz = sparse(1:n_sensors,1:n_sensors,nz_diag,n_sensors,n_sensors);

I1 = mu0/(4*pi)*(Ny*Rz-Nz*Ry);
I2 = mu0/(4*pi)*(Nz*Rx-Nx*Rz);
I3 = mu0/(4*pi)*(Nx*Ry-Ny*Rx);

end

