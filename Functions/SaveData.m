function cfg = SaveData(h, psiI, Q, fluxI, cfg, fm)
    % Function that saves the relevant information. This function saves under
    % the following convention
    %
    % The save name consists of several values. 
    % 
    % First number: trial number. Every seperate run has its own trial number
    %
    % Second number: save number. Each time it saves this is added to. We only
    % keep the 3 most recent copies where the rest are automatically deleted.
    %
    % Save String: The string inputted by the user to reference.
    %
    % R: River flux on/off
    % C: CSG flux on/off
    % P: Pumping on/off
    % Q: Evapotranspiration on/off
    % t: Years (rounded down) since simulation began. 
    %
    % Author: Samuel Dudley

    saveString = cfg.saveString;
    t_total = cfg.t_total;

    addpath(genpath('Data'))

    % Get a list of the files.
    files = dir(fullfile('Data', '*.mat'));

    % Get a list of the file names.
    fileCell = {files.name};

    % Key that seperates trial number to rest of name. 
    key = '_';

    % List to save all old save numbers in.
    savej = [];

    % If the trial number already exists we simply take it to be trialMax
    % instead and don't check in the loop.
    if isfield(cfg, 'trialNumber')
        trialNumber = cfg.trialNumber;
        % Don't need to check for more
        check1 = false;
    else
        check1 = true;
        trialMax = 0; % Initial guess.
    end

    if ~isfield(cfg, 'saveNumber') % Maybe it was deleted on accident?
        % We look to find the saveNumber if it does not exist here. 
        check2 = true;
    else
        saveNumber = cfg.saveNumber + 1; % add one to the save number
        check2 = false;
    end

    % fprintf('\ntrialNumber = %g\n', trialNumber)

    % Loop over all the files, check whether the trial number is maximum.
    if check1 
        for j = 1:length(fileCell)

            % Find the location of the key.
            keyLoc = strfind(fileCell{j}, key);
            
            % We need to make sure we only search over those that have data.
            if ~isempty(keyLoc)
                % Save the trial number
                trialj = str2double(fileCell{j}(1:keyLoc(1)-1));
                % Take the maximum between all previous and currently found. 
                trialMax = max(trialMax, trialj(1));
            end
        end
        % Set current trial
        cfg.trialNumber = trialMax + 1;
        trialNumber = cfg.trialNumber;
    end

    % Get a list of all save names
    for j = 1:length(fileCell)
        keyLoc = strfind(fileCell{j}, key);
        if ~isempty(keyLoc)
            if str2double(fileCell{j}(1:keyLoc(1)-1)) == trialNumber
                savej = [savej, str2double(fileCell{j}(keyLoc(1)+1:keyLoc(2)-1))];
            end
        end
    end


    % Set the current saveNumber
    if check2
        saveNumber = max([savej, 0]) + 1;
    end
    % Save new saveNumber.
    cfg.saveNumber = saveNumber;

    % Value for savename.
    t_total = floor(t_total/365); % convert to years (round down)

    % String save name.
    saveName = sprintf('%g_%g_%s_R%g_C%g_P%g_Q%g_t%.2fyears.mat', ...
        trialNumber, ...
        saveNumber, ...
        saveString, ...
        cfg.riverFlux, ...
        cfg.csgFlux, ...
        cfg.pumpOn, ...
        cfg.evapo, ...
        t_total);

    % This checks that a certain trial isn't being overwritten
    if any(savej > saveNumber)
        fprintf('\nWarning!\n\n You may be overwriting files.\nPlease check and then press anykey to continue...\n\n')
        pause
    end

    % Full save output.
    saveOutput = fullfile('Data', saveName);

    % Save data
    save(saveOutput, 'cfg', 'h', 'psiI', 'fluxI', 'Q', 'fm')

    % Say we have saved
    fprintf('\nSaved successfully at %g iterations\n\n', cfg.iter)


    % Now we delete old ones. 


    % we only want to delete those that are greater than 2 saves ago. 
    if isfield(cfg, 'saveDelete')
        if cfg.saveDelete
            toDelete = savej(savej+2 < saveNumber);
            if ~isempty(toDelete)
                fprintf('Deleting old saves\n')
            end
            
            % Search over all these, create the string FIX LATER
            for j = 1:length(toDelete)
                delString = sprintf('%g_%g_*.mat', trialNumber, toDelete(j));
                delString = fullfile('Data', delString);
                fileDelete = dir(delString);
                cellDelete = {fileDelete.name};
                if length(cellDelete) ~= 1
                    warning('Files may not have deleted correctly')
                end
                delete(fullfile('Data', cellDelete{1}))
            end
            
        end
    end

end





