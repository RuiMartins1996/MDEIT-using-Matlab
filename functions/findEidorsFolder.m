function eidorsFolder = findEidorsFolder(currentDirectory)

%Find EIDORS version
d = dir(currentDirectory);
isub = [d(:).isdir]; % returns logical vector

nameFolds = {d(isub).name}'; %cell array containing the folder names
nameFolds(ismember(nameFolds,{'.','..'})) = [];

% Search for folder containing the string eidors
str = '\w*eidors\w*';
matchStr = [];
for id = 1:length(nameFolds)
    matchStr = regexp(nameFolds{id},str);
    if not(isempty(matchStr))
        break
    end
end

if isempty(matchStr)
    error(strcat('Could not find EIDORS folder in ',currentDirectory));
end

eidorsFolder = nameFolds{id};

end