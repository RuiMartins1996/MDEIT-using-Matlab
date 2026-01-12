function modelParametersArray = ...
    sweep_model_parameters(parameters, matrix, maxValues, numSteps,model_parameters,mode)

if nargin <6
    mode = 'linear';
end

% If parameters is a char array (single string), convert to 1x1 cell array
if ischar(parameters)
    parameters = {parameters};
end

% Number of parameters
nParams = numel(parameters);

if nParams>1
    % Interpret minValues as a matrix with n_rows = nParams and n_cols
    % equal to the number of model_parameters to be generated
    parameter_values = matrix;

    modelParametersArray = ...
        sweep_model_parameters_multiple_parameter(parameters, parameter_values ,model_parameters);

else
    minValues = matrix;
    modelParametersArray = ...
        sweep_model_parameters_single_parameter(parameters, minValues, maxValues, numSteps,model_parameters,mode);
end


end




function     modelParametersArray = ...
    sweep_model_parameters_multiple_parameter(parameters, parameter_values ,model_parameters)

% Optional modelParameters template
if ~isempty(model_parameters)
    template = model_parameters;
    template = parseModelParameters(template);
else
    template = standardParameters();
end

% Number of parameters
nParams = numel(parameters);

assert(size(parameter_values,1) == nParams,'Expected number of rows in parameter_values to be nParams');

% Preallocate struct array
N = size(parameter_values,2);
modelParametersArray(1, N) = template;

% Loop to fill struct array
for i = 1:N
    temp = template;
    for k = 1:nParams
        temp.(parameters{k}) = parameter_values (k,i);
    end
    modelParametersArray(i) = temp;
end


end





function modelParametersArray = ...
    sweep_model_parameters_single_parameter(parameters, minValues, maxValues, numSteps,model_parameters,mode)

% Optional modelParameters template
if nargin > 4 && ~isempty(model_parameters)
    template = model_parameters;
    template = parseModelParameters(template);
else
    template = standardParameters();
end

% Number of parameters
nParams = numel(parameters);
assert(nParams == 1,'Expected number of parameters to be 1');

% Input checks
if ~(numel(minValues)==nParams && numel(maxValues)==nParams && numel(numSteps)==nParams)
    error('All input arrays must have the same length as parameters');
end

% Create vectors of sweep values for each parameter
sweepVectors = cell(1, nParams);
for k = 1:nParams

    if strcmp(mode,'linear')
        sweepVectors{k} = linspace(minValues(k), maxValues(k), numSteps(k));
    elseif strcmp(mode,'log')
        sweepVectors{k} = logspace(log10(minValues(k)), log10(maxValues(k)), numSteps(k));
    else
        error('mode not recognized');
    end

    % Check field exists in standardParameters
    if ~isfield(standardParameters(), parameters{k})
        error('Field "%s" does not exist in modelParameters.', parameters{k});
    end
end

% Generate all combinations using ndgrid
[grid{1:nParams}] = ndgrid(sweepVectors{:});

% Flatten each grid into a vector
for k = 1:nParams
    grid{k} = grid{k}(:);
end

% Total number of combinations
N = numel(grid{1});

% Preallocate struct array
modelParametersArray(1, N) = template;

% Loop to fill struct array
for i = 1:N
    temp = template;
    for k = 1:nParams
        temp.(parameters{k}) = grid{k}(i);
    end
    modelParametersArray(i) = temp;
end

end