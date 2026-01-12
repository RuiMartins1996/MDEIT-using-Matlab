function modelParameters = parseModelParameters(modelParameters)

[isValid, missingFields, extraFields] = checkModelParameters(modelParameters);

if ~isValid
    stdParameters = standardParameters();
    assert(~isempty(modelParameters), 'modelParameters cannot be empty.');

    % Handle missing fields
    if ~isempty(missingFields)
        for k = 1:numel(missingFields)
            modelParameters.(missingFields{k}) = stdParameters.(missingFields{k});
        end
    end

    % Handle extra fields
    if ~isempty(extraFields)
        warning('Removing unexpected field(s) from modelParameters: %s', strjoin(extraFields, ', '));
        modelParameters = rmfield(modelParameters, extraFields);
    end

    % Reorder to match standardParameters
    modelParameters = orderfields(modelParameters, stdParameters);
end

assert(strcmp(modelParameters.stimulationType,'adjacent'),'stimulationType must be adjacent');

end

function [isValid, missingFields, extraFields] = checkModelParameters(modelParameters)
%CHECKMODELPARAMETERS Validate modelParameters against standardParameters.
%
%   [isValid, missingFields, extraFields] = checkModelParameters(modelParameters)
%
%   Returns:
%     isValid        - logical true if all required fields exist and no extras
%     missingFields  - cell array of missing field names
%     extraFields    - cell array of unexpected field names

    template = standardParameters();

    expectedFields = fieldnames(template);
    providedFields = fieldnames(modelParameters);

    % Identify discrepancies
    missingFields = setdiff(expectedFields, providedFields);
    extraFields   = setdiff(providedFields, expectedFields);

    % Model parameters are valid only if both are empty
    isValid = isempty(missingFields) && isempty(extraFields);
end