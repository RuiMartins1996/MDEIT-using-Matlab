function [model_parameters_array_new,fmdl_cell_array] = ...
    mk_mdeit_model(model_parameters_array,model_folder,options)

% Parse options
%TODO
if nargin<3
    options = [];
end

opts = parse_options(options);

if nargin<2
    model_folder = [];
else
    if not(isempty(model_folder))
        assert(ischar(model_folder) || isstring(model_folder))
    end
end

%%

assert( (isstruct(model_parameters_array) && isvector(model_parameters_array)) || ...
    (iscell(model_parameters_array) && isvector(model_parameters_array)), ...
    'modelParametersCell must be a 1D struct array or a 1D cell array.');

N = length(model_parameters_array);

% Convert modelParametersArray to struct array
if iscell(model_parameters_array)
    model_parameters_array = cell2mat(model_parameters_array(:));  % forces 1D column first
    model_parameters_array = reshape(model_parameters_array, 1, []); % make 1×N row
end
%% Create multiple forward models

model_parameters_array_new = repmat(standardParameters(), 1, N);

% Preallocate 1xN struct array
fmdl_cell_array = cell(1, N);

for i = 1:length(model_parameters_array)
    model_parameters = model_parameters_array(i);

    fprintf('Building/Loading model %i of %i\n',i,N)

    [model_parameters,fmdl] = ...
        mk_mdeit_model_single_modelParameters(model_parameters,model_folder,opts);

    model_parameters_array_new(i) = model_parameters;
    fmdl_cell_array{i} = fmdl;
end

clc;

end


%% Function: Create a single fmdl from modelParameters
function [model_parameters,fmdl] = ...
    mk_mdeit_model_single_modelParameters(model_parameters,model_folder,opts)
%% Parse model parameters to correct missing fields
model_parameters = parseModelParameters(model_parameters);

%% Create file name with hashing for storing in disk
model_name = model_parameters_to_file_name(model_parameters,model_folder); %empty if modelFolder is empty

%% Distinguish between mesh and fmdl. Some fmdls have the same mesh, but different model_parameters
relevant_fields = {'maxsz','isCylindrical','height','radius','material','is2D','electrodeRadius','numOfRings','numOfElectrodesPerRing'};

% Initialize empty struct
mesh_parameters = struct();

% Copy only relevant fields from model_parameters
for k = 1:numel(relevant_fields)
    field = relevant_fields{k};
    if isfield(model_parameters, field)
        mesh_parameters.(field) = model_parameters.(field);
    end
end

mesh_name = model_parameters_to_mesh_name(mesh_parameters,model_folder);

%% If the mesh does not exist, generate it and save it
if not(exist(mesh_name, 'file') == 2) || opts.recompute == true
    %% Create EIDORS forward model
    if isfield(model_parameters,'material') && numel(fieldnames(model_parameters.material))>0 % Has material
        fmdl = create_model_with_material(model_parameters);
    else % No material
        fmdl= create_model(model_parameters);
    end

    %% Save into file
    if ~isempty(mesh_name)
        save(mesh_name,'fmdl');
    end
else % If this file already exists, load it
    var = load(mesh_name);
    fmdl = var.fmdl;
    mesh_parameters = file_name_to_mesh_parameters(mesh_name,model_folder);
end

%% If this model does not exist, generate it and save
if not(exist(model_name, 'file') == 2) || opts.recompute == true
    %% Assign the electrode contact impedances (since the EIDORS algorithms generate them with default values)
    fmdl = assign_contact_impedances(fmdl,model_parameters);
    
    %% Create stimulation patterns
    switch model_parameters.stimulationType
        case 'adjacent'
            stim = mk_stim_patterns(model_parameters.numOfElectrodesPerRing,model_parameters.numOfRings,[0,1],[0,1],{'meas_current'},1);
            fmdl.stimulation = stim;
        otherwise
            error(strcat('Stimulation type',model_parameters.stimulationType,' has not been implemented!'));
    end

    %% Place magnetometers
    if opts.no_magnetometers == true
        fprintf('Skipping magnetometer assignment \n');
    else
        fprintf('Assign magnetometers \n');
        fmdl = assign_magnetometers(fmdl,model_parameters);
    end
    %% Compute geometry matrices
    mu0 = model_parameters.mu0;

    if opts.no_geometry_matrices == true
        fprintf('Skipping computation of geometry matrices \n');
    else
        fprintf('Computing geometry matrices\n');
        fmdl = compute_geometry_matrices(fmdl,mu0);
    end

    %img.fwd_solve.get_all_meas = 1;

    %% Save the file name in a fmdl field
    fmdl.file_name = model_name;

    %% Save into file
    if ~isempty(model_name)
        save(model_name,'fmdl');
    end

    %% If this file already exists, load it
else
    var = load(model_name);
    fmdl = var.fmdl;
    model_parameters = file_name_to_model_parameters(model_name,model_folder);
end

end


%% Function: create_model_with_material
function fmdl = create_model_with_material(model_parameters)
    % The problem of netgen hanging is usually due to the orthobrick
    % having exact intersection with the domain (intersecting a boundary of
    % the domain)
    if model_parameters.is2D == true
        fmdl = create_model_with_material_2d(model_parameters);
    else
        fmdl = create_model_with_material_3d(model_parameters);
    end
end
%% Function: create_model_with_material_3d
function fmdl = create_model_with_material_3d(model_parameters)

dh = model_parameters.height/(model_parameters.numOfRings+1);
ring_vert_pos = dh:dh:model_parameters.height-dh;

electrode_radius = model_parameters.electrodeRadius;
maxsz = model_parameters.maxsz;

position = model_parameters.material.position;
radius = model_parameters.material.radius;
name = model_parameters.material.name;
type = model_parameters.material.type;

cyl_shape = [model_parameters.height,model_parameters.radius,model_parameters.maxsz];
elec_pos = [model_parameters.numOfElectrodesPerRing,ring_vert_pos];
% elec_shape = [electrode_radius, 0, maxsz ]; %circular electrodes
elec_shape = [electrode_radius, electrode_radius, maxsz ]; %square electrodes


if strcmp(type,'cylindrical')
    % Cylindrical anomaly (in constructive solid geometry language (CSG)
    % the cylinder is infinite, so it must be cut by intersection with
    % orthobrick)
    string1 = sprintf('cylinder(%.2f,%.2f,%.2f;%.2f,%.2f,%.2f;%.2f)',...
        position(1),position(2),-5*model_parameters.height,...
        position(1),position(2),5*model_parameters.height,...
        radius);

    % Bounding box of the cylindrical domain
    string2 = sprintf('orthobrick(%.2f,%.2f,%.2f;%.2f,%.2f,%.2f)',...
        -model_parameters.radius,-model_parameters.radius,0,...
        model_parameters.radius,model_parameters.radius,model_parameters.height);

    string_total = sprintf('solid ball = %s and %s -maxh=%.2f;',string1,string2,maxsz);

    % I don't think the problem has to do with maxsz in
    % string_total, since when it is constant from mesh to mesh,
    % the problem still only occurs in some meshes. Furthermore,
    % not setting -maxh still gives an error, so the problem is
    % elsewhere.
    % string_total = sprintf('solid ball = %s and %s;',string1,string2);

    % For some reason, it seems to work when the cylinder is fully
    % contained in the domain. That is, the bounding box that cuts
    % the cylinder is a bit shorther than the domain

    extra={'ball',string_total};

    fmdl= ng_mk_cyl_models(cyl_shape,elec_pos,elec_shape,extra);
    
elseif strcmp(type,'spherical')
    
    % Sphere in CSG
    string = sprintf('sphere(%.2f,%.2f,%.2f;%.2f)',...
        position(1),position(2),position(3),radius);

    string_total = sprintf('solid ball = %s -maxh=%.2f;',string,maxsz);
    
    extra={'ball',string_total};
    
    fmdl= ng_mk_cyl_models(cyl_shape,elec_pos,elec_shape,extra);
else
    error('type "%s" does not make sense or is not implemented',type)
end

%Sanity check to see if the material was inserted
assert(numel(fmdl.mat_idx) == 2,'fmdl.mat_idx should have length 2');

end
%% Function: create_model_with_material_2d
function fmdl = create_model_with_material_2d(model_parameters)

maxsz = model_parameters.maxsz;

position = model_parameters.material.position;
radius = model_parameters.material.radius;
name = model_parameters.material.name;
type = model_parameters.material.type;

electrode_radius = model_parameters.electrodeRadius;

cyl_shape = [0,model_parameters.radius,model_parameters.maxsz];
elec_shape = [electrode_radius, 0, maxsz ];

if strcmp(type,'cylindrical')
    tolerance = model_parameters.radius/100;
    % Cylindrical anomaly (in constructive solid geometry language (CSG)
    % the cylinder is infinite, so it must be cut by intersection with
    % orthobrick)
    string1 = sprintf('cylinder(%.2f,%.2f,%.2f;%.2f,%.2f,%.2f;%.2f)',...
        position(1),position(2),-5*model_parameters.height,...
        position(1),position(2),5*model_parameters.height,...
        radius);
    
    % Bounding box of the cylindrical domain
    string2 = sprintf('orthobrick(%.2f,%.2f,%.2f;%.2f,%.2f,%.2f)',...
        -model_parameters.radius,-model_parameters.radius,0,...
        model_parameters.radius,model_parameters.radius,1);

    string_total = sprintf('solid ball = %s and %s -maxh=%.2f;',string1,string2,maxsz);

    extra={'ball',string_total};

    fmdl= ng_mk_cyl_models(cyl_shape,model_parameters.numOfElectrodesPerRing,elec_shape,extra);

    %Sanity check to see if the material was inserted
    assert(numel(fmdl.mat_idx) == 2,'fmdl.mat_idx should have length 2');
else
    error('type "%s" does not make sense or is not implemented',type)
end

end
%% Function: create_model
function fmdl = create_model(model_parameters)
if model_parameters.is2D == true
    fmdl = create_model_2d(model_parameters);
else
    fmdl = create_model_3d(model_parameters);
end
end
%% Function: create_model_3d
function fmdl = create_model_3d(model_parameters)

dh = model_parameters.height/(model_parameters.numOfRings+1);
ring_vert_pos = dh:dh:model_parameters.height-dh;
electrode_radius = model_parameters.electrodeRadius;
maxsz = model_parameters.maxsz;

cyl_shape = [model_parameters.height,model_parameters.radius,model_parameters.maxsz];
elec_pos = [model_parameters.numOfElectrodesPerRing,ring_vert_pos];
elec_shape = [electrode_radius, electrode_radius, maxsz ];

fmdl= ng_mk_cyl_models(cyl_shape,elec_pos,elec_shape);

end

%% Function: create_model_2d
function fmdl = create_model_2d(model_parameters)

maxsz = model_parameters.maxsz;

electrode_radius = model_parameters.electrodeRadius;

cyl_shape = [0,model_parameters.radius,model_parameters.maxsz];
elec_shape = [electrode_radius, 0, maxsz ];

fmdl= ng_mk_cyl_models(cyl_shape,model_parameters.numOfElectrodesPerRing,elec_shape);
end

%% Function: 

function fmdl = assign_contact_impedances(fmdl,model_parameters)
    
    assert(isfield(fmdl,'electrode'),'Expected an .electrode field in fmdl')
    
    for i = 1:numel(fmdl.electrode)
        fmdl.electrode(i).z_contact = model_parameters.electrodeContactImpedance;
    end

end
%% Function: assignMagnetometers
function fmdl = assign_magnetometers(fmdl,model_parameters)

num_sensors = model_parameters.numOfSensors;
sensor_radius = model_parameters.sensorRadius;
measurement_axis_type = model_parameters.measurementAxisType;

dim = size(fmdl.nodes,2);

% Override sensor placement function if the sensor positions have been
% specified
if isfield(model_parameters,'sensorPositions')
    if not(isempty(model_parameters.sensorPositions))
        sensor_positions = model_parameters.sensorPositions;
        if dim == 2 || dim == 3
            assert(size(sensor_positions,1) == num_sensors);

            sensor_locations =  sensor_positions;

            sensor_axes = ...
                repmat(struct('axis1', [], 'axis2', [],'axis3',[]), 1, num_sensors);

            for m = 1:num_sensors
                sensor_axes(m).axis1 = [1,0,0];
                sensor_axes(m).axis2 = [0,1,0];
                sensor_axes(m).axis3 = [0,0,1];
            end

            % Assign sensorLocations and sensorAxes to fmdl
            sensors = repmat(struct('position', [], 'axes', []), 1, num_sensors);

            for m = 1:num_sensors
                sensors(m).position = sensor_locations(m,:);
                sensors(m).axes = sensor_axes(m);
            end

            fmdl.sensors = sensors;
            
            return
        end
    end
end


% TODO: There should be a model_parameters parameter called configuration.
% It can be 'spherical', 'cylindrical' and 'plane' or others. If the model
% is 2D, then 'spherical and 'cylindrical' are incompatible and there
% should be an error. But for now lets do it like this for backwards
% compatibility

if dim == 2
    options{1} = num_sensors; %number of magnetometers per ring is simply the number of sensors
    options{2} = 1; %number of rings is 1 in 2D
    options{3} = 0; %the height of the sensors should be 0
    options{4} = sensor_radius; %the radius is correct

    [sensor_locations,sensor_axes]  =...
        place_magnetometers(options,'cylindrical',measurement_axis_type);
else
    switch model_parameters.isCylindrical
        case true

            options{2} = model_parameters.numOfRings;

            numMagnetPerRing = num_sensors/model_parameters.numOfRings;

            isInteger = abs(numMagnetPerRing - round(numMagnetPerRing)) < 1e-10;

            if not(isInteger)
                error('Number of sensors must be a multiple of number of sensors rings')
            end

            options{1} = numMagnetPerRing;

            dh = model_parameters.height/(model_parameters.numOfRings+1);
            heights = dh:dh:model_parameters.height-dh;

            options{3} = heights;
            options{4} = sensor_radius;

            [sensor_locations,sensor_axes]  =...
                place_magnetometers(options,'cylindrical',measurement_axis_type);
        case false
            options{1} = num_sensors;
            options{2} = sensor_radius;
            options{3} = sum(fmdl.nodes,1)/length(fmdl.nodes);

            [sensor_locations,sensor_axes]  =...
                place_magnetometers(options,'spherical',measurement_axis_type);
    end

end

% Assign sensorLocations and sensorAxes to fmdl
sensors = repmat(struct('position', [], 'axes', []), 1, num_sensors);

for m = 1:num_sensors
    sensors(m).position = sensor_locations(m,:);
    sensors(m).axes = sensor_axes(m);
end

fmdl.sensors = sensors;
end

%% FUNCTION: Parse options
function opts = parse_options(options)
%PARSE_OPTIONS Convert a cell array of option strings to a struct of flags.
%
% Usage:
%   opts = parse_options(options, valid_options)
%
% Inputs:
%   options        - cell array of option strings, e.g. {'no_geometry_matrices','no_magnetometers'}
% Output:
%   opts           - struct with logical fields for each valid option

valid_options = {'no_geometry_matrices','no_magnetometers','recompute'};

% Default: all flags off
opts = struct();
for i = 1:numel(valid_options)
    opts.(valid_options{i}) = false;
end

% If no options provided, just return defaults
if nargin < 1 || isempty(options)
    return;
end

% Validate and enable specified options
for i = 1:numel(options)
    opt = lower(strtrim(options{i}));
    if ismember(opt, valid_options)
        opts.(opt) = true;
    else
        warning('Unknown option: "%s"', opt);
    end
end
end