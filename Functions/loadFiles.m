%% Load files.
files = dir(fullfile('Data', ['*' loadString '*.mat']));

fileCell = {files.name};

if length(fileCell) < 1
    error('The file string you are searching for does not currently exist')
end

if length(fileCell) > 1
    fprintf('\nMultiple files found.\n')
    fprintf('Please select from the following\n\n')
    for j = 1:length(fileCell)
        fprintf('\tID %g: %s\n', j, fileCell{j})
    end
    prompt = '\nPlease select a file ID Number:  ';
    fileID = input(prompt);

    while fileID > length(fileCell)
        prompt = 'Please select a valid file ID Number:  ';
        fileID = input(prompt);
    end

    fileToLoad = fileCell{fileID};
    fprintf('\n\nYou chose: %s\n\n', fileToLoad)
else
    fileToLoad = fileCell{1};
    fprintf('\n\nLoading: %s\n\n', fileToLoad)
end

if strcmpi(loadString, 'FullSS_UniformMesh.mat')
    warning('Make sure  set a new saveString variable in cfgnew')
end

clear loadString prompt j fileID files fileCell



load(fullfile('Data', fileToLoad))