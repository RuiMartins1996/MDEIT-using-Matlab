function modelParameters = ...
    createEITandMDEITCylindricalModel(modelFolder,modelParameters)

folderName = modelFolder;

% Check if modelParameters has every expected field
try
    checkModelParameters(modelParameters);
    
    % Unpack the valid structure
    fields = fieldnames(modelParameters);
    for k = 1:numel(fields)
        % create a variable with the same name as the field
        % in THIS function workspace
        eval([fields{k} ' = modelParameters.(fields{k});']);
    end

catch %if model parameters does not have every field

    stdParameters = standardParameters();

    % Copy existing fields
    assert(not(isempty(modelParameters)));

    fields = fieldnames(modelParameters);
    for k = 1:numel(fields)
        field = fields{k};
        stdParameters.(field) = modelParameters.(field);
    end
    
    modelParameters = stdParameters;
end

assert(strcmp(modelParameters.stimulationType,'adjacent'),'stimulationType must be adjacent');

%% Create file name with hashing
[hashStr, shortHash] = encodeModelParameters(modelParameters);

% Store mapping from shortHas to long hasStr
mappingFile = fullfile(modelFolder, 'hashMapping.mat');
if exist(mappingFile, 'file')
    load(mappingFile, 'hashMapping');
else
    hashMapping = containers.Map;
end
hashMapping(shortHash) = hashStr;
save(mappingFile, 'hashMapping');


fileName = strcat(...
    strcat(folderName,'/model_',shortHash,'.mat'));

%% Sanity check
S = decodeModelParameters(hashStr);

if not(compareStructs(S, modelParameters))
    error('Decoded struct and modelParameters are not the same!')
end

%% If this file does not exist, generate it and save
if not(exist(fileName, 'file') == 2)

    dh = modelParameters.height/(modelParameters.numOfRings+1);
    ring_vert_pos = dh:dh:modelParameters.height-dh;
    
    rng(modelParameters.randomConductivitySeed);

    %% Create  forward model
    fmdl= ng_mk_cyl_models(...
        [modelParameters.height,modelParameters.radius,modelParameters.maxsz],[modelParameters.numOfElectrodesPerRing,ring_vert_pos],[0.4,0.4,modelParameters.maxsz]);
    
    %% Create stimulation patterns
    switch modelParameters.stimulationType
        case 'adjacent'
            stim = mk_stim_patterns(modelParameters.numOfElectrodesPerRing,modelParameters.numOfRings,[0,1],[0,1],{'meas_current'},1);
            fmdl.stimulation = stim;
        otherwise
            error(strcat('Stimulation type',modelParameters.stimulationType,' has not been implemented!'));
    end

    %% Create homogeneous model
    img = mk_image(fmdl,1.0);
    %% Create random conductivity distribution if required
    if modelParameters.hasRandomConductivity

        rng(modelParameters.randomConductivitySeed);
        
        minCond = 0.9;
        maxCond = 1.1;

        % Random conductivities between a and b
        img.elem_data = minCond + (maxCond-minCond).*rand(size(img.elem_data,1),1);
    end

    %% Create anomaly if required
    if modelParameters.hasRandomAnomaly
        %true if modelParameters has field 'anomaly', override LEGACY behaviour
        if isfield(modelParameters,'anomaly')
            img = mk_image(fmdl,1.0);
            img = placeAnomaly(img,modelParameters);
        else
            % Place anomaly if requested (LEGACY)
            img = mk_image(fmdl,1.0);
            img = placeAnomaly(img);
        end
    end

    %% Place magnetometers
   
    numOfSensors = modelParameters.numOfSensors;
    
    R = modelParameters.sensorRadius;

    if not(modelParameters.isCylindrical)           % Place magnetometers spherically
        
        options{1} = numOfSensors;
        options{2} = R;
        options{3} = sum(fmdl.nodes,1)/length(fmdl.nodes);
        sensorLocations = placeMagnetometers(options,'spherical');
    else                            % Place magnetometers cilindrically
        
        numOfSensorRings = modelParameters.numOfRings;
        options{2} = numOfSensorRings;

        numMagnetPerRing = numOfSensors/numOfSensorRings;
        
        isInteger = abs(numMagnetPerRing - round(numMagnetPerRing)) < 1e-10;
        
        if not(isInteger)
            error('Here')
        end

        options{1} = numMagnetPerRing;

        heights = linspace(0,modelParameters.height,numOfSensorRings);

        options{3} = heights;

        options{4} = R;

        sensorLocations = placeMagnetometers(options,'cylindrical');
    end

    img.sensorLocations = sensorLocations;

    img.fwd_solve.get_all_meas = 1;
    %% Compute geometry matrices
    mu0 = modelParameters.mu0;
    img = computeGeometryMatrices(img,mu0);

    %% Save into file
    save(fileName,'img');
end
end

%% FUNCTION: placeSphericalAnomaly
function img = placeSphericalAnomaly(img,anomalyConductivity,anomalyRadius,anomalyCenter)

baseConductivity = 1.0;

str = strcat('(x-',num2str(anomalyCenter(1)),...
    ').^2 + (y-',num2str(anomalyCenter(2)),...
    ').^2 + (z-',num2str(anomalyCenter(3)),...
    ').^2 <',num2str(anomalyRadius),'^2');

select_fcn = inline(str,'x','y','z');

memb_frac = elem_select(img.fwd_model, select_fcn);
img = mk_image(img.fwd_model, baseConductivity + (anomalyConductivity-baseConductivity)*memb_frac );
end

%% FUNCTION: findSphericalAnomalyCenter

function anomalyCenter = findSphericalAnomalyCenter(modelParameters)
% Find an admisseable center to place an anomaly of anomalyRadius

anomalyRadius = modelParameters.anomaly.radius;

numOfRings = modelParameters.numOfRings;

dh = modelParameters.height/(numOfRings+1);
ring_vert_pos = dh:dh:modelParameters.height-dh;

%Anomaly must be inside imaging plane, so ring_vert_pos must be
%accounted!
upperBound = max(ring_vert_pos);
lowerBound = min(ring_vert_pos);

%Ensure that it will fit inside the imaging plane
ub = upperBound-anomalyRadius;
lb = lowerBound+anomalyRadius;

if ub<lb
    error('Anomaly will not fit in this model')
end

%Generate random cylindrical coordinates 
absAnomalyZ = lb + (ub - lb) * rand(1);

ub = modelParameters.radius-anomalyRadius;
lb = 0;

if ub<lb
    error('Anomaly will not fit in this model')
end
absAnomalyR = lb + (ub - lb) * rand(1);

absAnomalyTheta = 2*pi*rand(1); 

% Convert to cartesian coordinates
absAnomalyX = absAnomalyR*cos(absAnomalyTheta);
absAnomalyY = absAnomalyR*sin(absAnomalyTheta);

absAnomalyCenter = [absAnomalyX;absAnomalyY;absAnomalyZ];

sign = randi([0,1],3,1) * 2 - 1; % Produces -1 or +1

% Compute anomaly center
anomalyCenter = sign.*absAnomalyCenter;
end
%% FUNCTION: placeAnomaly
function img = placeAnomaly(img, varargin)

if nargin == 1 %Default use
    anomalyConductivity = 3.0;
    anomalyRadius = img.modelParameters.radius/5;

    anomalyCenter = findSphericalAnomalyCenter(img,anomalyRadius);
    img = placeSphericalAnomaly(img,anomalyConductivity,anomalyRadius,anomalyCenter);
else %If struct modelParameters is given
    modelParameters = varargin{1};
    
    if not(isfield(modelParameters,'anomaly'))
        error('No field "anomaly" found for second argument');
    end

    switch modelParameters.anomaly.type
        case 'spherical'
            anomalyCenter = findSphericalAnomalyCenter(modelParameters);
            img = placeSphericalAnomaly(img,modelParameters.anomaly.conductivity,modelParameters.anomaly.radius,anomalyCenter);
        otherwise
            error(strcat('Type',string(modelParameters.anomaly.type),' is not implemented'));
    end
end

end



%% FUNCTION: checkModelParameters
% Check if the structure modelParameters is admisseable to be parsed by the
%main function
function checkModelParameters(modelParameters)
expectedFields = {'maxsz', 'isCylindrical', 'height','radius', ...
    'hasRandomAnomaly','anomaly','hasRandomConductivity','randomConductivitySeed', 'numOfRings', ...
    'numOfElectrodesPerRing', 'stimulationType','numOfSensors','sensorRadius','mu0'};

for k = 1:numel(expectedFields)
    if ~isfield(modelParameters, expectedFields{k})
        error('Missing required field: %s', expectedFields{k});
    end
end



end
%% FUNCTION: standardParameters
function stdParameters = standardParameters()

stdParameters = struct();

stdParameters.maxsz = 0.5;
stdParameters.isCylindrical = true;
stdParameters.height = 3;
stdParameters.radius = 1;

stdParameters.hasRandomAnomaly = false;
stdParameters.anomaly = struct();
stdParameters.hasRandomConductivity = false;
stdParameters.randomConductivitySeed = 0;

stdParameters.numOfRings = 2;
stdParameters.numOfElectrodesPerRing = 4;
stdParameters.stimulationType = 'adjacent';
stdParameters.mu0 = 1;

% Assume same number of measurements for every stimulation pattern!
stdParameters.numOfSensors = 16;

stdParameters.sensorRadius = 2.0;

end

function isSame = compareStructs(a, b)
    % Compare two structs recursively, ignoring row/column orientation
    if ~isstruct(a) || ~isstruct(b)
        error('Inputs must be structs');
    end

    % Compare field names
    if ~isequal(sort(fieldnames(a)), sort(fieldnames(b)))
        isSame = false;
        return;
    end

    % Loop over fields
    flds = fieldnames(a);
    for i = 1:numel(flds)
        f = flds{i};
        va = a.(f);
        vb = b.(f);

        if isstruct(va) && isstruct(vb)
            % Recursively compare structs
            if ~compareStructs(va, vb)
                isSame = false;
                return;
            end
        elseif isnumeric(va) && isnumeric(vb)
            % Compare numerics, ignoring row/column shape
            if ~isequal(va(:), vb(:))
                isSame = false;
                return;
            end
        elseif ischar(va) && ischar(vb)
            if ~strcmp(va, vb)
                isSame = false;
                return;
            end
        elseif isstring(va) && isstring(vb)
            if ~strcmp(va, vb)
                isSame = false;
                return;
            end
        elseif islogical(va) && islogical(vb)
            if ~isequal(va(:), vb(:))
                isSame = false;
                return;
            end
        else
            % Fallback strict check
            if ~isequal(va, vb)
                isSame = false;
                return;
            end
        end
    end

    isSame = true;
end