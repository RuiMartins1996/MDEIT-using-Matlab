clc; clear all; close all;

%% Prepare workspace
% Get the full path of the current script
fullpath = mfilename('fullpath');
% Extract just the folder
script_folder = fileparts(fullpath);
cd(script_folder);

% Have to add the functions path manually so prepare_workspace runs
grandparent_folder =fileparts(fileparts(script_folder));
addpath(genpath(fullfile(grandparent_folder,'functions')));

model_folder = prepare_workspace(script_folder);

%% Setup EIDORS
clc;

rng(1)

data_folder = strcat(script_folder ,'\data');

%% Define the characteristic scales in SI units

z0 = 0.0058; %(Ohm m^2) is the contact impedance from the CEM article 58 Ohm cm^2
l0 = 40e-3; %(m) the tank radius
I0 = 2.4e-3;%(A) the magnitude of the injected current

% The derived characteristic units
V0 = z0*I0/(l0^2); %(V)
sigma0 = l0/z0; %(S/m)
J0 = I0/(l0^2);
H0 = I0/l0; 
B0 = 1.25663706127e-6*I0/l0; 

%% Build/load multiple forward models of different sensor radius

model_parameters = create_default_3d_model_parameters(l0, z0, sigma0, I0);

model_parameters.maxsz = model_parameters.radius/3;
model_parameters.material = struct();
model_parameters.numOfElectrodesPerRing = 8;
model_parameters.numOfRings = 4;

% model_parameters.numOfSensors = model_parameters.numOfRings*model_parameters.numOfElectrodesPerRing;

model_parameters.numOfSensors = model_parameters.numOfRings*12;

num_sensors = model_parameters.numOfSensors;
num_sensors_per_ring = num_sensors/model_parameters.numOfRings;
sensor_radius = model_parameters.radius;

dh = model_parameters.height/(model_parameters.numOfRings+1);
heights = dh:dh:model_parameters.height-dh;

%% Sweep parameters
r_min = 1.05*model_parameters.radius;
r_max = 3*model_parameters.radius;
n_steps = 20;

all_sensor_radius = linspace(r_min,r_max,n_steps);

s_circular.sensor_radius = all_sensor_radius;
s_square.sensor_radius = all_sensor_radius;
s_triangle.sensor_radius = all_sensor_radius;

%% Create cylindrical sensor positions

fmdl_circular_array = cell(n_steps,1);

for k = 1:n_steps
    sensor_radius = all_sensor_radius(k);
    
    for n = 1:numel(heights)
        height = heights(n);
        theta = 0:2*pi/num_sensors_per_ring:2*pi/num_sensors_per_ring*(num_sensors_per_ring-1);
        sensor_positions((n-1)*num_sensors_per_ring+1:n*num_sensors_per_ring,:) = [sensor_radius*cos(theta)',sensor_radius*sin(theta)',ones(num_sensors_per_ring,1)*height];
    end

    model_parameters.sensorPositions = sensor_positions;

    % Build models
    [model_parameters,fmdl_array] = mk_mdeit_model(model_parameters,model_folder,[]);
    fmdl_circular = fmdl_array{1};
    
    fmdl_circular_array{k} = fmdl_circular;

end



s_circular.model_parameters = model_parameters;
%% Create square sensor positions

fmdl_square_array = cell(n_steps,1);

for k = 1:n_steps

    sensor_radius = all_sensor_radius(k);

    % Preallocate
    sensor_positions = zeros(num_sensors, 3);

    % Square parameters
    L = 2 * sensor_radius;           % side length
    P = 4 * L;                       % perimeter
    s = linspace(0, P, num_sensors_per_ring + 1);
    s(end) = [];                     % remove duplicate endpoint

    for n = 1:numel(heights)
        height = heights(n);

        x = zeros(num_sensors_per_ring, 1);
        y = zeros(num_sensors_per_ring, 1);

        for i = 1:num_sensors_per_ring
            si = s(i);

            if si < L
                % Bottom edge (left → right)
                x(i) = -sensor_radius + si;
                y(i) = -sensor_radius;

            elseif si < 2*L
                % Right edge (bottom → top)
                x(i) =  sensor_radius;
                y(i) = -sensor_radius + (si - L);

            elseif si < 3*L
                % Top edge (right → left)
                x(i) =  sensor_radius - (si - 2*L);
                y(i) =  sensor_radius;

            else
                % Left edge (top → bottom)
                x(i) = -sensor_radius;
                y(i) =  sensor_radius - (si - 3*L);
            end
        end

        idx = (n-1)*num_sensors_per_ring + 1 : n*num_sensors_per_ring;
        sensor_positions(idx, :) = [x, y, height*ones(num_sensors_per_ring,1)];
    end

    model_parameters.sensorPositions = sensor_positions;

    % Build models
    [model_parameters,fmdl_array] = mk_mdeit_model(model_parameters,model_folder,[]);
    fmdl_square = fmdl_array{1};

    fmdl_square_array{k} = fmdl_square;

end


s_square.model_parameters = model_parameters;

%% Create triangular sensor positions


fmdl_triangle_array = cell(n_steps,1);

for k = 1:n_steps

    sensor_radius = all_sensor_radius(k);

    sensor_positions = zeros(num_sensors, 3);

    % === Triangle geometry (incircle radius = sensor_radius) ===
    r = sensor_radius;
    a = 2*sqrt(3)*r;

    V1 = [ -a/2, -r ];
    V2 = [  a/2, -r ];
    V3 = [   0 , 2*r ];

    vertices = [V1; V2; V3];

    for n = 1:numel(heights)
        height = heights(n);

        pts = vertices;   % start with vertices

        % Each edge is a cell array of segments
        edges = {
            [V1; V2]
            [V2; V3]
            [V3; V1]
            };

        % Refinement loop
        while size(pts,1) < num_sensors_per_ring
            new_edges = {};

            for e = 1:numel(edges)
                segments = edges{e};

                for s = 1:size(segments,1)/2
                    P = segments(2*s-1,:);
                    Q = segments(2*s,:);

                    R = 0.5*(P + Q);

                    % Add midpoint if we still need sensors
                    if size(pts,1) < num_sensors_per_ring
                        pts = [pts; R];
                    else
                        break;
                    end

                    % Subdivide segment
                    new_edges{end+1} = [P; R];
                    new_edges{end+1} = [R; Q];
                end

                if size(pts,1) >= num_sensors_per_ring
                    break;
                end
            end

            edges = new_edges;

            % Safety break (in case someone asks for too few sensors)
            if isempty(edges)
                break;
            end
        end

        % Trim in case we slightly overshot
        pts = pts(1:num_sensors_per_ring, :);

        idx = (n-1)*num_sensors_per_ring + 1 : n*num_sensors_per_ring;
        sensor_positions(idx,:) = [pts, height*ones(num_sensors_per_ring,1)];
    end

    model_parameters.sensorPositions = sensor_positions;

    % Build models
    [model_parameters,fmdl_array] = mk_mdeit_model(model_parameters,model_folder,[]);
    fmdl_triangle = fmdl_array{1};

    fmdl_triangle_array{k} = fmdl_triangle;

end



s_triangle.model_parameters = model_parameters;

%% Show figures side by side

figure
subplot(1,3,1)
hold on
show_fem(fmdl_circular);
plot_sensors(fmdl_circular);
hold off
grid on;grid minor;box on;
axis([-2*sensor_radius 2*sensor_radius -2*sensor_radius 2*sensor_radius])
view(2)
subplot(1,3,2)
hold on
show_fem(fmdl_square);
plot_sensors(fmdl_square);
axis([-2*sensor_radius 2*sensor_radius -2*sensor_radius 2*sensor_radius])
hold off
grid on;grid minor;box on;
view(2)
subplot(1,3,3)
hold on
show_fem(fmdl_triangle);
plot_sensors(fmdl_triangle);
axis([-2*sensor_radius 2*sensor_radius -2*sensor_radius 2*sensor_radius])
hold off
grid on;grid minor;box on;
view(2)
%% Compute the jacobian for each of them

condition_number_array_circular = nan(n_steps,3);
condition_number_array_square = nan(n_steps,3);
condition_number_array_triangle = nan(n_steps,3);

for i = 1:n_steps
    for d = 1:3
        select_sensor_axis = d;

        % Circular configuration
        imgh = mk_image_mdeit(fmdl_circular_array{i},1.0);
        lambdatimesdAdp = @(lambda) []; %legacy
        A = @(sigma) M(imgh,sigma);

        J = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',select_sensor_axis);

        singular_values = svds(J,rank(J));
        condition_number_array_circular(i,d) = singular_values(1)/singular_values(end);

        s_circular.file_name = create_file_name(data_folder,s_circular.model_parameters,'circular');
        s_circular.condition_number_array = condition_number_array_circular;
        save(s_circular.file_name,"s_circular");

        % Square configuration
        imgh = mk_image_mdeit(fmdl_square_array{i},1.0);
        lambdatimesdAdp = @(lambda) []; %legacy
        A = @(sigma) M(imgh,sigma);

        J = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',select_sensor_axis);

        singular_values = svds(J,rank(J));
        condition_number_array_square(i,d) = singular_values(1)/singular_values(end);

        s_square.file_name = create_file_name(data_folder,s_square.model_parameters,'square');
        s_square.condition_number_array = condition_number_array_square;
        save(s_square.file_name,"s_square");

        % Triangle configuration
        imgh = mk_image_mdeit(fmdl_triangle_array{i},1.0);
        lambdatimesdAdp = @(lambda) []; %legacy
        A = @(sigma) M(imgh,sigma);

        J = calc_jacobian_mdeit(imgh,imgh.elem_data,lambdatimesdAdp,A,'mdeit1',select_sensor_axis);

        singular_values = svds(J,rank(J));
        condition_number_array_triangle(i,d) = singular_values(1)/singular_values(end);

        s_triangle.file_name = create_file_name(data_folder,s_triangle.model_parameters,'triangle');
        s_triangle.condition_number_array = condition_number_array_triangle;
        save(s_triangle.file_name,"s_triangle");
    end

end


%% FUNCTIONS
function file_name = create_file_name(data_folder,model_parameters,string)
    name = sprintf('data_%s_E_%i_R_%i_M_%i',...
        string,...
        model_parameters.numOfElectrodesPerRing,...
        model_parameters.numOfRings,...
        model_parameters.numOfSensors);
    file_name = strcat(data_folder,"\",name);
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