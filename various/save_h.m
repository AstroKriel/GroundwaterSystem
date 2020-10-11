function save_h(cfg, h, t_total)

    addpath(genpath('Data'))

    files = dir(fullfile('Data', '*.mat'));

    fileCell = {files.name};

    key = '_';

    if isfield(cfg, 'trialNumber')
        check1 = false; 
        trialNumber = cfg.trialNumber;
    else
        check1 = true;
        trialMax = 0;
    end

    if check1 
        for j = 1:length(fileCell)
            
            keyLoc = strfind(fileCell{j}, key);
            
            if ~isempty(keyLoc)
                % Save the trial number
                trialj = str2double(fileCell{j}(1:keyLoc(1)-1));
                % Take the max between all previous found ones, and this. 
                trialMax = max(trialMax, trialj(1));
            end
        end
        % Set the current trial into the cfg.
        cfg.trialNumber = trialMax + 1; 
        trialNumber = cfg.trialNumber;
    end

    save_hString = fullfile('Analysis', [num2str(trialNumber), '_', cfg.saveString, '.mat']);

    if isfile(fullfile(save_hString))
        % File exists
        matFileObj = matfile(save_hString);
        matFileObj.Properties.Writable = true;
        
        matFileObj.h = [matFileObj.h, h];
        matFileObj.t_total = [matFileObj.t_total, t_total];
    else
        save(save_hString, 'h', 't_total', 'cfg');
    end
    
end
