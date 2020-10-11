function cfg = loadcfg(newcfg, cfg)
    % Function combines new and old cfg.
    
    if isfield(newcfg, 'Solvek')
        if strcmpi(newcfg.Solvek, 'upwinding')
            cfg.Solvek = str2func('Upwinding');
        elseif strcmpi(newcfg.Solvek, 'harmonic')
            cfg.Solvek = str2func('Harmonick');
        elseif strcmpi(newcfg.Solvek, 'arithmetic')
            cfg.Solvek = str2func('Arithmetic');
        else
            error('Solvek incorrect')
        end
    end

    if isfield(newcfg, 'Jacobian')
        if strcmpi(newcfg.Jacobian, 'Banded')
            cfg.Jacobian = str2func('JacobianBanded');
        elseif strcmpi(newcfg.Jacobian, 'Inexact')
            cfg.Jacobian = str2func('JacobianInexact');
        else
            error('Jacobian incorrect')
        end
    end

    if isfield(newcfg, 'rainfallModel')
        if strcmpi(newcfg.rainfallModel, 'constant')
            cfg.rainFunc = str2func('rain_constant');
        elseif strcmpi(newcfg.rainfallModel, 'simple')
            cfg.rainFunc = str2func('rain_basic'); 
        elseif strcmpi(newcfg.rainfallModel, 'complex')
            cfg.rainFunc = str2func('rain_random'); 
        else
            error('Rainfall model is incorrect')
        end
    end

    if isfield(newcfg, 'yearRain')
        cfg.yearRain = newcfg.yearRain;
    end

    % Drought information
    if isfield(newcfg, 'droughtLength')
        % drought start is required
        if ~isfield(newcfg, 'droughtStart')
            error('For drought length, a start year is required')
        elseif newcfg.droughtStart*365 <= cfg.t_total
            error('droughtStart is less than current time.')
        end
        
        cfg.droughtStart = newcfg.droughtStart*365;
        cfg.droughtLength = newcfg.droughtLength*365;
    elseif isfield(newcfg, 'droughtStart')
        % Does not necessarily need a length
        if newcfg.droughtStart*365 <= cfg.t_total
            error('droughtStart is less than current time.')
        else
            cfg.droughtStart = newcfg.droughtStart*365;
            cfg.droughtLength = inf; % set to infinity.
        end
    end
    
    if isfield(newcfg, 't_total')
        cfg.t_total = newcfg.t_total;
    end

    if isfield(newcfg, 'dt_max')
        cfg.dt_max = newcfg.dt_max;
    end

    if isfield(newcfg, 'dt_min')
        cfg.dt_min = newcfg.dt_min;
    end

    if ((isfield(newcfg, 'dt_max') && newcfg.dt > cfg.dt_max) ...
            || (isfield(newcfg, 'dt_min') && newcfg.dt < cfg.dt_min)) ...
            && isfield(newcfg, 'dt')
        error('initial time step must be within bounds set')
    elseif isfield(newcfg, 'dt')
        % set initial time step
        cfg.dt = newcfg.dt;
    end

    if isfield(newcfg, 't_max')
        if newcfg.t_max * 365 - 10 < cfg.t_total
            error('t_max is too close to, or less than current time, please increase it')
        else
            % Convert from years to days (easier). 
            cfg.t_max = 365 * newcfg.t_max;
        end
    end
    
    if isfield(newcfg, 'debugmode')
        error('cannot change debugmode')
    end
    
    % Efficient 
    if isfield(newcfg, 'efficient')
        cfg.efficient = newcfg.efficient;
    end

    % Required saveString
    if isfield(newcfg, 'saveString')
        cfg.saveString = newcfg.saveString;
    end

    % At what iteration modular do we save results?
    if isfield(newcfg, 'saveITER')
        cfg.saveITER = newcfg.saveITER;
    end

    if isfield(newcfg, 'theta_1')
        if newcfg.theta_1 ~= 1 && newcfg.theta_1 ~= 1/2
            error('theta_1 and theta_2 must be either 1 or 0.5')
        else
            cfg.theta_1 = newcfg.theta_1;
        end
    end

    if isfield(newcfg, 'theta_2')
        if newcfg.theta_2 ~= 1 && newcfg.theta_2 ~= 1/2
            error('theta_1 and theta_2 must be either 1 or 0.5')
        else
            cfg.theta_2 = newcfg.theta_2;
        end
    end
    
    if isfield(newcfg, 'riverFlux')
        cfg.riverFlux = newcfg.riverFlux;
    end
    if isfield(newcfg, 'csgFlux')
        cfg.csgFlux = newcfg.csgFlux;
    end
    if isfield(newcfg, 'rainFlux')
        cfg.rainFlux = newcfg.rainFlux;
    end
    if isfield(newcfg, 'pumpOn')
        cfg.pumpOn = newcfg.pumpOn;
    end
    if isfield(newcfg, 'evapo')
        cfg.evapo = newcfg.evapo;
    end
    
    if ~isfield(cfg, 'pumpX')
        [~, pumpX] = min(abs(cfg.x - 100));
        pumpZl = find(cfg.z >= 75, 1, 'last');
        pumpZu = find(cfg.z >= 55, 1, 'last');
        cfg.pumpNl = (pumpX - 1) * cfg.r + pumpZl; % node 1
        cfg.pumpNu = (pumpX - 1) * cfg.r + pumpZu; % node 2
    end
    
    % Set up so these turn on or off at given times.
    if isfield(newcfg, 'csgFluxTime')
        cfg.csgFluxTime = newcfg.csgFluxTime .* 365;
    end

    if isfield(newcfg, 'rainFluxTime')
        cfg.rainFluxTime = newcfg.rainFluxTime .* 365;
    end

    if isfield(newcfg, 'pumpOnTime')
        cfg.pumpOnTime = newcfg.pumpOnTime .* 365;
    end
    
    if isfield(newcfg, 'K_R')
        cfg.K_R = newcfg.K_R;
    end 
    
    if isfield(newcfg, 'x_CSG')
        cfg.x_CSG = newcfg.x_CSG;
    end
    
    if isfield(cfg, 'trialNumber')
        cfg = rmfield(cfg, 'trialNumber');
    end
    
    if isfield(cfg, 'saveNumber')
        cfg = rmfield(cfg, 'saveNumber');
    end
    
    if isfield(newcfg, 'showPlots')
        cfg.showPlots = newcfg.showPlots;
    end
    
    if isfield(newcfg, 'plotSizing')
        cfg.plotSizing = newcfg.plotSizing;
    end
    
    if isfield(newcfg, 'saveH')
        cfg.saveH = newcfg.saveH;
    end
    
    if isfield(newcfg, 'saveDelete')
        cfg.saveDelete = newcfg.saveDelete;
end