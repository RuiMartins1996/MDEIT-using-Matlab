function model_parameters = file_name_to_model_parameters(file_name,model_folder)
% Extract shortHash from fileName
tokens = regexp(file_name, 'model_([^.]+)\.mat', 'tokens');
shortHash = tokens{1}{1};

% Get shortHash to longHash mapping
mappingFile = fullfile(model_folder, 'hash_mapping.mat');
if exist(mappingFile, 'file')
    load(mappingFile, 'hash_mapping');
else
    error('hash_mapping should exist')
end

hash_str = hash_mapping(shortHash);

model_parameters = decodeModelParameters(hash_str);

end

