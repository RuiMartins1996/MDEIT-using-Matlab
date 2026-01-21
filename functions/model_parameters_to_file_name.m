function file_name = model_parameters_to_file_name(model_parameters,model_folder)

% Return an empty file name if modelFolder was not specified
if isempty(model_folder)
    file_name = [];
    return;
end

[hash_str, short_hash_str] = encodeModelParameters(model_parameters);

% Store mapping from short_hash to long hash_str
mapping_file = fullfile(model_folder, 'hash_mapping.mat');
if exist(mapping_file, 'file')
    load(mapping_file, 'hash_mapping');
else
    hash_mapping = containers.Map;
end
hash_mapping(short_hash_str) = hash_str;
save(mapping_file, 'hash_mapping');

file_name = strcat(...
    strcat(model_folder,'\model_',short_hash_str,'.mat'));

% Sanity check
S = decodeModelParameters(hash_str);

if not(compare_structs(S, model_parameters))
    error('Decoded struct and modelParameters are not the same!')
end
end


%% FUNCTION: compareStructs
function is_same = compare_structs(a, b)
% Compare two structs recursively, ignoring row/column orientation
if ~isstruct(a) || ~isstruct(b)
    error('Inputs must be structs');
end

% Compare field names
if ~isequal(sort(fieldnames(a)), sort(fieldnames(b)))
    is_same = false;
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
        if ~compare_structs(va, vb)
            is_same = false;
            return;
        end
    elseif isnumeric(va) && isnumeric(vb)
        % Compare numerics, ignoring row/column shape
        if ~isequal(va(:), vb(:))
            is_same = false;
            return;
        end
    elseif ischar(va) && ischar(vb)
        if ~strcmp(va, vb)
            is_same = false;
            return;
        end
    elseif isstring(va) && isstring(vb)
        if ~strcmp(va, vb)
            is_same = false;
            return;
        end
    elseif islogical(va) && islogical(vb)
        if ~isequal(va(:), vb(:))
            is_same = false;
            return;
        end
    else
        % Fallback strict check
        if ~isequal(va, vb)
            is_same = false;
            return;
        end
    end
end

is_same = true;
end