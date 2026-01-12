%% FUNCTION: standardParameters
function stdParameters = standardParameters()

stdParameters = struct();

% Mesh resolution
stdParameters.maxsz = 0.5;

% Domain geometry
stdParameters.isCylindrical = false;
stdParameters.measurementAxisType = 'cartesian';
stdParameters.height = 3;
stdParameters.radius = 1;

stdParameters.material = struct();
stdParameters.is2D = false;

stdParameters.hasRandomAnomaly = false;
stdParameters.anomaly = struct();
stdParameters.hasRandomConductivity = false;
stdParameters.randomConductivitySeed = 0;

% Electrode configuration and properties
stdParameters.electrodeRadius = 0.4;
stdParameters.electrodeContactImpedance = 1;
stdParameters.numOfRings = 2;
stdParameters.numOfElectrodesPerRing = 4;
stdParameters.stimulationType = 'adjacent';

% Magnetic sensor configuration and magnetic permeability
stdParameters.mu0 = 1.2566370612720e-6 ;
stdParameters.numOfSensors = 16;
stdParameters.sensorRadius = 2.0;
stdParameters.sensorPositions = [];

end