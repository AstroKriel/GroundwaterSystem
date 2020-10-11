function [fm, cfg] = checkCFG(in, cfgMesh)
    % Takes inputs and constructs the cfg file computing all relevant
    % calculations needed. 

    % Initialise cfg file:
    cfg = load('FixedParameters.mat');

    if isfield(in, 'Solvek')
        if strcmpi(in.Solvek, 'upwinding')
            cfg.Solvek = str2func('Upwinding');
        elseif strcmpi(in.Solvek, 'harmonic')
            cfg.Solvek = str2func('Harmonick');
        elseif strcmpi(in.Solvek, 'arithmetic')
            cfg.Solvek = str2func('Arithmetic');
        else
            error('Solvek incorrect')
        end
    else
        % Set to default
        cfg.Solvek = str2func('Upwinding');
    end

    if isfield(in, 'Jacobian')
        if strcmpi(in.Jacobian, 'Banded')
            cfg.Jacobian = str2func('JacobianBanded');
        elseif strcmpi(in.Jacobian, 'Inexact')
            cfg.Jacobian = str2func('JacobianInexact');
        else
            error('Jacobian incorrect')
        end
    else
        % Set to default
        cfg.Jacobian = str2func('JacobianBanded'); 
    end

    if isfield(in, 'rainfallModel')
        if strcmpi(in.rainfallModel, 'constant')
            cfg.rainFunc = str2func('rain_constant');
        elseif strcmpi(in.rainfallModel, 'simple')
            fprintf('Cosine model chosen')
            pause(3)
            cfg.rainFunc = str2func('rain_basic'); 
        elseif strcmpi(in.rainfallModel, 'complex')
            cfg.rainFunc = str2func('rain_random'); 
        else
            error('Rainfall model is incorrect')
        end
    else
        % set to default
        cfg.rainFunc = str2func('rain_basic');
    end
    
    if isfield(in, 'yearRainAverage')
        cfg.yearRainAverage = in.yearRainAverage;
    else
        error('Yearly rain needed')
    end

    % Drought information
    if isfield(in, 'droughtLength')
        % drought start is required
        if ~isfield(in, 'droughtStart')
            error('For drought length, a start year is required')
        end
        
        cfg.droughtStart = in.droughtStart*365;
        cfg.droughtLength = in.droughtLength*365;
    elseif isfield(in, 'droughtStart')
        % Does not necessarily need a length
        cfg.droughtStart = in.droughtStart*365;
        cfg.droughtLength = inf; % set to infinity.
    end

    if ~isfield(in, 'dt_max')
        % default
        cfg.dt_max = 10;
    else
        cfg.dt_max = in.dt_max;
    end

    if ~isfield(in, 'dt_min')
        % default
        cfg.dt_min = 1 / (24 * 60);
    else
        cfg.dt_min = in.dt_min;
    end

    if (isfield(in, 'dt_max') && in.dt > cfg.dt_max) || (isfield(in, 'dt_min') && in.dt < cfg.dt_min)
        error('initial time step must be within bounds set')
    else
        % set initial time step
        cfg.dt = in.dt;
    end

    if isfield(in, 't_max')
        % Convert from years to days (easier). 
        cfg.t_max = 365 * in.t_max;
    else
        cfg.t_max = inf;
    end


    % Debugging mode
    if isfield(in, 'debugMode')
        cfg.debugMode = in.debugMode;
    else
        cfg.debugMode = false;
    end

    % Efficient 
    if isfield(in, 'efficient')
        cfg.efficient = in.efficient;
    else
        cfg.efficient = true;
    end

    % Required saveString
    if isfield(in, 'saveString')
        cfg.saveString = in.saveString;
    else
        error('Save string required')
    end

    % At what iteration modular do we save results?
    if isfield(in, 'saveITER')
        cfg.saveITER = in.saveITER;
    else
        cfg.saveITER = 100;
    end

    if isfield(in, 'theta_1')
        if in.theta_1 ~= 1 && in.theta_1 ~= 1/2
            error('theta_1 and theta_2 must be either 1 or 0.5')
        else
            cfg.theta_1 = in.theta_1;
        end
    else
        cfg.theta_1 = 1;
    end

    if isfield(in, 'theta_2')
        if in.theta_2 ~= 1 && in.theta_2 ~= 1/2
            error('theta_1 and theta_2 must be either 1 or 0.5')
        else
            cfg.theta_2 = in.theta_2;
        end
    else
        cfg.theta_2 = 1;
    end

    % Check flux values all given
    if ~isfield(in, {'riverFlux', 'csgFlux', 'rainFlux', 'pumpOn', 'evapo'})
        % Mitch: is this correct? Or is it only throwing this error if ALL of them are missing,
        % and ignores the case where one is there and the others aren't for example.
        error('all flux boolean parameters are required')
    else
        cfg.riverFlux = in.riverFlux;
        cfg.csgFlux = in.csgFlux;
        cfg.rainFlux = in.rainFlux;
        cfg.pumpOn = in.pumpOn;
        cfg.evapo = in.evapo;
    end

    % Set up so these turn on or off at given times.
    if isfield(in, 'csgFluxTime')
        cfg.csgFluxTime = in.csgFluxTime .* 365;
    else
        % If no information provided it turns on and never turns off.
        if in.csgFlux
            cfg.csgFluxTime = [0 inf];
        else
            cfg.csgFluxTime = [inf inf];
        end
    end

    if isfield(in, 'rainFluxTime')
        cfg.rainFluxTime = in.rainFluxTime .* 365;
    else
        if in.rainFlux
            cfg.rainFluxTime = [0 inf];
        else
            cfg.rainFluxTime = [inf inf];
        end
    end

    if isfield(in, 'pumpOnTime')
        cfg.pumpOnTime = in.pumpOnTime .* 365;
    else
        if in.pumpOn
            cfg.pumpOnTime = [0 inf];
        else
            cfg.pumpOnTime = [inf inf];
        end
    end

    % Set up mesh
    [fm, cfg] = setupMesh(cfgMesh, cfg);

    % Unpack mesh results
    if cfg.debugMode
        cfg.NM(~isnan(cfg.NM)) = 1; % Set all materials to be type 1;
    end
    cfg.NM0 = cfg.NM;
    cfg.NM0(isnan(cfg.NM)) = 5; % change all NAN values to a 5.

    % Find node that the river is on
    cfg.Nriver = find(cfg.z == cfg.R_b);

    % Adjust river perm
    if isfield(in, 'K_R')
        cfg.K_R = in.K_R;
    end % no need for an else statement because default is already set.

    % Find node that the CSG plant is on.
    if isfield(in, 'x_CSG')
        % Mitch: implement in.x_CSG - should be quick, but want to double-check my logic here
        % before I type it in.
        cfg.NCSG = cfg.N - cfg.r + find(cfg.z == 5);
        cfg.x_CSG = in.x_CSG;
    else
        cfg.NCSG = cfg.N - cfg.r + find(cfg.z == 5);
        cfg.x_CSG = 5000; % Set to 5km away
    end

    % Find the node that the pump is on. 
    % Mitch: check this for non-uniform - not 100% correct atm
    % Will need to fix nodes at these points within the non-uniform mesh I think?

    % Mitch: pump can be anywhere (not on a node) - check below for non-uniform
    % double-check this implementation
    [~, pumpX] = min(abs(cfg.x - 100));
    pumpZl = find(cfg.z >= 75, 1, 'last');
    pumpZu = find(cfg.z >= 55, 1, 'last');
    cfg.pumpNl = (pumpX - 1) * cfg.r + pumpZl; % node 1
    cfg.pumpNu = (pumpX - 1) * cfg.r + pumpZu; % node 2

    % Initialise iteration
    cfg.iter = 1; 

    % Are we showing plots?
    if isfield(in, 'showPlots')
        cfg.showPlots = in.showPlots;
    else
        cfg.showPlots = true;
    end

    % Plot window sizing & positioning
    if isfield(in, 'plotSizing')
        cfg.plotSizing = in.plotSizing;
    else
        cfg.plotSizing = [0.1, 0.05, 0.8, 0.9];
    end
    
    if isfield(in, 'P_R')
        cfg.P_R = in.P_R;
    end
    
    if isfield(in, 'saveH')
        cfg.saveH = in.saveH;
    else
        cfg.saveH = true;
    end
    
    if isfield(in, 'saveDelete')
        cfg.saveDelete = in.saveDelete;
    else
        cfg.saveDelete = false;
    end

end







