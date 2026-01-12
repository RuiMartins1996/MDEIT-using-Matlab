function models_folder = prepare_workspace(script_folder)
%PREPARE_WORKSPACE  Set up paths, folders, and global parameters.
%
%   prepare_workspace(script_folder)
%
%   INPUT:
%       script_folder — absolute path of the folder where the calling script lives
%
%   This function:
%       • Creates and adds a ./data folder inside script_folder
%       • Moves two directories up from script_folder
%       • Adds all subfolders of "functions" and "libraries"
%       • Runs globalParameters.m
%       • Creates and adds ./models folder
%

%% Move to script folder
if nargin < 1 || ~isfolder(script_folder)
    error('prepare_workspace:InvalidFolder', ...
        'Input script_folder must be a valid folder path.');
end
cd(script_folder);

%Create + add data folder
data_folder = fullfile(script_folder, 'data');
if ~exist(data_folder, 'dir')
    mkdir(data_folder);
end
addpath(data_folder);

%Add functions folder
functions_folder = find_in_parent_folders(script_folder, 'functions', 'folder');
if ~isempty(functions_folder)
    addpath(genpath(functions_folder));
else
    error('Could not find "functions" folder in parent directories.');
end

%Add libraries folder
libraries_folder = find_in_parent_folders(script_folder, 'libraries', 'folder');
if ~isempty(libraries_folder)
    addpath(genpath(libraries_folder));
else
    error('Could not find "libraries" folder in parent directories.');
end

% Set or create model folder
models_folder = find_in_parent_folders(script_folder, 'models', 'folder');
if ~isempty(models_folder)
    addpath(genpath("models"));
else
    error('Could not find "models" folder in parent directories.');
end

% Run globalParameters.m
gp_file = find_in_parent_folders(script_folder, 'globalParameters.m', 'file');
if ~isempty(gp_file)
    run(gp_file);
end

% Setup EIDORS
setupEidors(script_folder);

clc;


end

function found_path = find_in_parent_folders(start_folder, target_name, type)
%FIND_UPWARDS   Search upward through parent directories for a file or folder.
%
%   found_path = find_upwards(start_folder, target_name, type)
%
%   INPUTS:
%       start_folder : initial directory where the search begins
%       target_name  : name of file or folder to search for
%       type         : 'file' or 'folder'
%
%   OUTPUT:
%       found_path   : absolute path to the matched item, or '' if not found
%
%   The function searches parent directories up to 10 levels above.
%

% Validate
if nargin < 3
    error('find_upwards:NotEnoughInputs', ...
        'Usage: find_upwards(start_folder, target_name, type).');
end
if ~isfolder(start_folder)
    error('find_upwards:InvalidStart', 'start_folder must exist.');
end
if ~ischar(target_name) && ~isstring(target_name)
    error('find_upwards:InvalidName', ...
        'target_name must be a string or char.');
end
if ~any(strcmp(type, {'file','folder'}))
    error('find_upwards:InvalidType', ...
        'type must be ''file'' or ''folder''.');
end

% Initialization
current = start_folder;
max_levels = 10;
target_name = char(target_name);  % convert to char if needed

for lvl = 1:max_levels

    switch type
        case 'folder'
            candidate = fullfile(current, target_name);
            if isfolder(candidate)
                found_path = candidate;
                return;
            end

        case 'file'
            candidate = fullfile(current, target_name);
            if exist(candidate, 'file')
                found_path = candidate;
                return;
            end
    end

    % Step one directory up
    parent = fileparts(current);

    % Root reached
    if strcmp(parent, current)
        break;
    end

    current = parent;
end

% Not found
found_path = '';
end
