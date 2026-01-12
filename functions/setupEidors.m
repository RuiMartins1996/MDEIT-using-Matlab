function eidorsFolder = setupEidors(start_folder)
%SETUPEIDORS  Finds the EIDORS folder by searching upward and runs startup.m
%
%   eidorsFolder = setupEidors(start_folder)
%
%   Searches upward through parent directories starting from start_folder
%   until it finds a folder whose name matches "*eidors*".
%
%   Then:
%       • Adds EIDORS to the path
%       • Removes overloads path
%       • Runs eidors/startup.m
%
%   Returns:
%       eidorsFolder — full absolute path to the detected EIDORS folder.
%

    %% Validate input
    if nargin < 1 || ~isfolder(start_folder)
        error('setupEidors:InvalidInput', ...
              'start_folder must be a valid directory.');
    end

    %% ==== Search upward for an EIDORS directory ====
    eidorsFolder = find_eidors_folder(start_folder);

    if isempty(eidorsFolder)
        error('setupEidors:NotFound', ...
            'Could not locate an EIDORS installation in parent directories.');
    end

    fprintf('Found EIDORS at: %s\n', eidorsFolder);

    %% ==== Add EIDORS to path ====
    addpath(genpath(eidorsFolder));

    %% ==== Remove overloads (important for absolute solvers) ====
    overload_path = fullfile(eidorsFolder, 'eidors', 'overloads');
    if isfolder(overload_path)
        rmpath(genpath(overload_path));
    end

    %% ==== Run startup.m ====
    startupFile = fullfile(eidorsFolder, 'eidors', 'startup.m');
    if exist(startupFile, 'file')
        run(startupFile);
    else
        warning('setupEidors:StartupMissing', ...
            'startup.m was not found at %s', startupFile);
    end
end


%% ------------------------------------------------------------------------
function eidorsFolder = find_eidors_folder(start_folder)
% helper: searches upward for a folder containing "eidors" in its name

    max_levels = 10;
    current = start_folder;

    for k = 1:max_levels

        % List subfolders
        d = dir(current);
        isub = [d.isdir];
        subdirs = {d(isub).name};
        subdirs = subdirs(~ismember(subdirs,{'.','..'}));

        % Look for "*eidors*" (case-insensitive)
        idx = find(contains(lower(subdirs), 'eidors'));

        if ~isempty(idx)
            eidorsFolder = fullfile(current, subdirs{idx(1)});
            return;
        end

        % Move one level up
        parent = fileparts(current);
        if strcmp(parent, current)
            break; % reached filesystem root
        end
        current = parent;
    end

    % Nothing found
    eidorsFolder = '';
end




% function eidorsFolder = setupEidors(mainPath)
% 
% %setupEidors: Finds the EIDORS folder and calls the startup file. This
% %function assumes that mainPath is in the same directory as the EIDORS
% %folder.
% 
% %Change directory
% cd(mainPath);
% 
% %Find EIDORS version
% d = dir(cd);
% isub = [d(:).isdir]; % indexes of directories
% 
% nameFolds = {d(isub).name}'; %cell array containing the folder names
% nameFolds(ismember(nameFolds,{'.','..'})) = []; %remove . and .. dirs
% 
% % Search for folder containing the string eidors
% str = '\w*eidors\w*';
% for id = 1:length(nameFolds)
%     matchStr = regexp(nameFolds{id},str);
%     if not(isempty(matchStr))
%         break
%     end
% end
% 
% eidorsFolder = nameFolds{id};
% 
% addpath(genpath(eidorsFolder));
% rmpath(genpath(strcat(eidorsFolder,"/eidors/overloads"))) %THIS IS CAUSING AN ERROR IN THE ABSOLUTE SOLVER IF NOT REMOVED
% 
% startupFile = strcat(eidorsFolder,"/eidors/startup.m");
% run(startupFile)
% 
% end
% 
