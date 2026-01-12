function model_parameters = create_kai_3d_model_parameters(l0, z0, sigma0, I0)
model_parameters.maxsz = 5e-3/l0;

model_parameters.is2D = false;

model_parameters.isCylindrical = true;
model_parameters.height = 70e-3/l0;
model_parameters.radius = 40e-3/l0;

model_parameters.numOfRings = 1;
model_parameters.numOfElectrodesPerRing = 16;

model_parameters.electrodeRadius = 8e-3/l0;
model_parameters.electrodeContactImpedance = 0.0058/z0;

model_parameters.randomConductivitySeed = 0;

model_parameters.numOfSensors = 16;
model_parameters.sensorRadius = 45e-3/l0; 
model_parameters.mu0 = 1;

anomaly_conductivity = 1e-12/sigma0;
anomaly_position     = [20,0,35]*1e-3/l0;
anomaly_radius       = 25/2*1e-3/l0;

model_parameters.anomaly = struct( ...
    'type', 'spherical', ...
    'conductivity', anomaly_conductivity, ...
    'radius', anomaly_radius, ...
    'position', anomaly_position);
end
