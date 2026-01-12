function model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0)

%% Define forward model (2D real tank experiment)

model_parameters.maxsz = 5e-3/l0;

model_parameters.is2D = false;

model_parameters.isCylindrical = true;
model_parameters.height = 120e-3/l0;
model_parameters.radius = 40e-3/l0;

model_parameters.numOfRings = 1;
model_parameters.numOfElectrodesPerRing = 16;

model_parameters.electrodeRadius = 10e-3/l0; %10e-3/l0 é o que costuma ser
model_parameters.electrodeContactImpedance = 0.0058/z0;

model_parameters.numOfSensors = 16;% Changed here to 22 because of realistic opm noise case
model_parameters.sensorRadius = 70e-3/l0; % Changed to 45 to increase the signal strength
model_parameters.mu0 = 1;

anomaly_conductivity = 1e-12/sigma0;
anomaly_position     = [20,0,60]*1e-3/l0;
anomaly_radius       = 25/2*1e-3/l0;

model_parameters.anomaly = struct( ...
    'type', 'cylindrical', ...
    'conductivity', anomaly_conductivity, ...
    'radius', anomaly_radius, ...
    'position', anomaly_position);

model_parameters.material = struct( ...
    'type', 'cylindrical', ...
    'name', 'plastic_cylinder', ...
    'radius', anomaly_radius, ...
    'position', anomaly_position);

end
