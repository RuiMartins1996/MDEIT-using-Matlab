function [sensor_locations,sensor_axes]  = place_magnetometers(options,configuration_type,measurement_axis_type)

if nargin<2
    configuration_type = 'cylindrical';
end

if nargin<3
    measurement_axis_type = 'cartesian';
end

valid_configuration_types = {'cylindrical','spherical'};
assert(...
    ismember(configuration_type,valid_configuration_types),...
    sprintf('configuration_type: %s is not valid',configuration_type))

valid_measurement_axis_type = {'cartesian','cylindrical','spherical'};
assert(...
    ismember(measurement_axis_type,valid_measurement_axis_type ),...
    sprintf('measurement_axis_type: %s is not valid',measurement_axis_type))

opts = parse_options(options,configuration_type);

sensor_locations = place_magnetometers_positions(opts,configuration_type);

sensor_axes = construct_sensor_axis(sensor_locations,measurement_axis_type);

end

%% FUNCTIONS: parse_options
function opts = parse_options(options,configuration_type)
    switch configuration_type
        case 'cylindrical'
            opts = struct();
            opts.num_sensors_per_ring = options{1};
            opts.num_rings = options{2};
            opts.heights = options{3};
            opts.sensor_radius = options{4};
        case 'spherical'
            opts = struct();
            opts.num_sensors = options{1};
            opts.sensor_radius = options{2};
            opts.center = options{3};
    end
end

%% FUNCTIONS: place_magnetometers
function sensor_locations = place_magnetometers_positions(opts,configuration_type)

switch configuration_type
    case 'cylindrical'
        sensor_locations  = ...
            place_magnetometers_cylindrical( ...
            opts.num_rings,opts.num_sensors_per_ring,...
            opts.heights,opts.sensor_radius);

    case 'spherical'
        sensor_locations  = ...
            place_magnetometers_spherical(opts.num_sensors,opts.center,opts.sensor_radius);
end

end
%% FUNCTIONS: place_magnetometers_cylindrical
function sensor_locations  = place_magnetometers_cylindrical(...
    num_rings,...
    num_sensors_per_ring,...
    heights,...
    sensor_radius)

num_sensors = num_sensors_per_ring*num_rings;
sensor_locations = zeros(num_sensors,3);

for m = 1:num_rings
    for n = 1:num_sensors_per_ring

        theta = n*2*pi/num_sensors_per_ring;
        sensor_id = n+(m-1)*num_sensors_per_ring;

        sensor_locations(sensor_id,1) = sensor_radius*cos(theta);
        sensor_locations(sensor_id,2) = sensor_radius*sin(theta);
        sensor_locations(sensor_id,3) = heights(m);

    end
end

end

%% FUNCTIONS: place_magnetometers_spherical
function sensor_locations  = place_magnetometers_spherical(...
    num_sensors,...
    center,...
    sensor_radius)

% Sample uniform points in the surface of the sphere
[V,~]=SpiralSampleSphere(num_sensors);

sensor_locations = sensor_radius*V+center;

end

%% FUNCTIONS: construct_sensor_axis
function sensor_axes = construct_sensor_axis(sensor_locations,measurement_axis_type)

num_sensors = size(sensor_locations,1);

center_of_referential = mean(sensor_locations);

sensor_axes = ...
    repmat(struct('axis1', [], 'axis2', [],'axis3',[]), 1, num_sensors);

switch measurement_axis_type
    case 'cartesian'
        for m=1:num_sensors
            sensor_axes(m).axis1 = [1,0,0];
            sensor_axes(m).axis2 = [0,1,0];
            sensor_axes(m).axis3 = [0,0,1];
        end
    case 'cylindrical'
        for m=1:num_sensors
            p = sensor_locations(m,:)-center_of_referential;

            x = p(1);
            y = p(2);
            z = p(3);

            rhat = [x,y,0];
            rhat = rhat/norm(rhat);

            thetahat = [-y,x,0];
            thetahat = thetahat/norm(thetahat);

            zhat = [0,0,1];


            sensor_axes(m).axis1 = rhat;
            sensor_axes(m).axis2 = thetahat;
            sensor_axes(m).axis3 = zhat;
        end
    case 'spherical'
        for m = 1:num_sensors
            
            p = sensor_locations(m,:)-center_of_referential;

            x = p(1);
            y = p(2);
            z = p(3);
            
            r = sqrt(x^2+y^2+z^2);
            theta = acos(z/sqrt(r^2));
            phi = sign(y)*acos(x/sqrt(x^2+y^2));

            rhat = sin(theta)*cos(phi)*[1,0,0]+sin(theta)*sin(phi)*[0,1,0]+cos(theta)*[0,0,1];
            thetahat = cos(theta)*cos(phi)*[1,0,0]+cos(theta)*sin(phi)*[0,1,0]-sin(theta)*[0,0,1];
            phihat = -sin(phi)*[1,0,0]+cos(phi)*[0,1,0];

            sensor_axes(m).axis1 = rhat;
            sensor_axes(m).axis2 = thetahat;
            sensor_axes(m).axis3 = phihat;
        end
   
    otherwise
        error('measurement_axis_type %s is not valid',measurement_axis_type);
end



end
