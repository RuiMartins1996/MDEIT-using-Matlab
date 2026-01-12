function img = createMDEITmodel(maxsz,numElectrodesPerRing,numRings,numOfSensors)
%% Global Parameters
baseConductivity = 1.0;
anomalyConductivity = 3.0;

height = 3;
radius = 1;

if nargin<2
    numElectrodesPerRing= 8;
    ring_vert_pos = [1/3,2/3,1,4/3,5/3,2,7/3,8/3];
    numRings = length(ring_vert_pos);

    numOfSensors = numElectrodesPerRing*numRings;
else
    dh = height/(numRings+1);
    ring_vert_pos = dh:dh:height-dh;
end

%% Create EIT forward model

% extra={'ball','solid ball = sphere(0,0,.5;.1);'};%sphere(x,y,z;radius)

fmdl= ng_mk_cyl_models(...
    [height,radius,maxsz],[numElectrodesPerRing,ring_vert_pos],[0.2,0.2,maxsz]);

% Stimulation patterns
stim = mk_stim_patterns(numElectrodesPerRing,numRings,[0,1],[0,1],{'meas_current'},1);
fmdl.stimulation = stim;


img = mk_image(fmdl,1.0);
numOfIterations = 0;

while all(abs(img.elem_data - 1)< 1e-12) && numOfIterations<10
    % Create image and add a circular anomaly at (rand,rand,rand) with radius 0.1
    anomalyRadius = radius/5;

    %Anomaly must be inside imaging plane, so ring_vert_pos must be
    %accounted!
    upperBound = max(ring_vert_pos);
    lowerBound = min(ring_vert_pos);

    ub = upperBound-2*anomalyRadius; %ensure they are a bit separate from the top and bottom of the cylinder too.
    lb = lowerBound+2*anomalyRadius;
    
    absAnomalyZ = lb + (ub - lb) * rand(1);
    
    ub = radius-2*anomalyRadius;
    lb = anomalyRadius/2;

    absAnomalyR = lb + (ub - lb) * rand(1);
    absAnomalyTheta = 2*pi*rand(1); %ensure (x,y) are points in a circle
    
    absAnomalyX = absAnomalyR*cos(absAnomalyTheta);
    absAnomalyY = absAnomalyR*sin(absAnomalyTheta);

    absAnomalyCenter = [absAnomalyX;absAnomalyY;absAnomalyZ];

    sign = randi([0,1],3,1) * 2 - 1; % Produces -1 or +1
    
    anomalyCenter = sign.*absAnomalyCenter;

    % anomalyCenter = num2str((1-2*anomalyRadius)*rand(3,1)+anomalyRadius);
    str = strcat('(x-',num2str(anomalyCenter(1,:)),...
        ').^2 + (y-',num2str(anomalyCenter(2,:)),...
        ').^2 + (z-',num2str(anomalyCenter(3,:)),...
        ').^2 <',num2str(anomalyRadius),'^2');

    select_fcn = inline(str,'x','y','z');

    memb_frac = elem_select(fmdl, select_fcn);
    img = mk_image(fmdl, baseConductivity + (anomalyConductivity-baseConductivity)*memb_frac );
    numOfIterations = numOfIterations +1;
end

if numOfIterations>=10
    error('Failed to create non-homogeneous model! This will result in errors in other algorithms!');
else
    fprintf('Succeeded in creating non-homgeneous model!');
    show_fem(img);
    pause(3)
end

img.fwd_solve.get_all_meas = 1;

%% Place magnetometers

R = 2.5;
options{1} = numOfSensors;
options{2} = R;
options{3} = sum(fmdl.nodes,1)/length(fmdl.nodes);
sensorLocations = placeMagnetometers(options,'spherical');

img.sensorLocations = sensorLocations;
end

