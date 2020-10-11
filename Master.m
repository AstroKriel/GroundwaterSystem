%% Master file for final implementation

% Formatting
clear; clc; close all
format short
format compact

% Add filepaths
addpath( genpath('Mesh') )
addpath( genpath('Functions') )
addpath( genpath('Rainfall') )
addpath( genpath('various') )

% Specify the string to load from.
loadString = 'SteadyState';

% From here on in,
% This file is the hub for all activity that occurs.

% First: Would you like to load from a previous save? 
load_OLD = false;

% Mitch: not 100% sure what is happening here will check soon
% I've set this to false for now, so not concerned atm.
if load_OLD
    % Now we construct the string to search for.
    loadFiles
end

% If you are loading from an old file, you can change the parameters
% to run with by defining new parameters in the newcfg struct.
newcfg = [];

newcfg.showPlots = 1;
newcfg.saveString = ''; % set this to something unique
newcfg.pumpOn = true;
newcfg.ksolve = 'Harmonic';

% For example, if you want to change any parameters (i.e. change whether flux is turned
% on or change theta_1 or theta_2 etc you can do so here. For example if you
% wanted to change theta_1 to crank-nicholson you would set:
% newcfg.theta_1 = 1/2;

% Analysis - from steady state
% newCFG stuff will be to load from there.

% The above changes will be combined with the old cfg here.
if load_OLD
    % Combine with newcfg
    cfg = loadcfg(newcfg, cfg);
end

%% If NEW RUN

% Is this a fresh run? Ie, we're not loading from an old file.
if ~load_OLD
    % We need to specify parameters for cfg.
    in = [];
    cfgMesh = [];
    
    % --- Flux Information --- %
    in.riverFlux = true;
    in.csgFlux = false;

    % Time in years that csgFluxTime is active (defaults to immediately)
    in.csgFluxTime = [5, 8];
    in.rainFlux = true;
    
    % --- Flow Rate Information --- %
    in.pumpOn = false;
    in.pumpOnTime = [5 inf]; % years
    in.evapo = true;

    % --- River river hydraulic conductivity (m/d) --- %
    in.K_R = 0.5; % 0.5 is default

    % --- CSG Location --- %
    % in.x_CSG = 5000; % Using default.
    
    % --- Solution methodology --- %
    % in.Jacobian = 'Inexact'; % Inexact or [Banded]
    % in.Solvek = 'Upwinding'; % [Upwinding] or Arithmetic or Harmonic
    
    % --- Rainfall model --- %
    % 'simple'      (cosine model)
    % 'constant'    (constant rainfall)
    % 'complex'     (least squares cosine model)
    
    in.yearRainAverage = 0.44; % metres/year
    in.rainfallModel = 'constant';
    
    % --- Droughts --- %
    in.droughtStart = 5; % starts at 5 years. 
    in.droughtLength = 5; % length in years.
    
    % --- Other parameters --- %
    in.dt     = 0.5; % Initial time step [days]
    in.dt_max = 15;  % Maximum time step [days]
    % in.dt_min = 1/(24*60); % Minimum time step (1 min is default).
    in.t_max  = inf; % t_max in years. Set to inf for indefinite. 

    in.theta_1 = 1; % theta value for flux calculations
    in.theta_2 = 1; % theta value for Flow Rate calculations
    
    % --- Saving parameters --- %
    in.saveString = 'MitchTesting';
    in.saveITER = 50; % save every ~ iterations. Default = 100.
    
    % --- Mesh parameters --- %
    cfgMesh.r = 21; % Base number of rows
    cfgMesh.c = 51; % Base number of columns
    cfgMesh.r_div = 1; % # of row subdivisions
    cfgMesh.c_div = 1; % # of column subdivisions
    cfgMesh.type = 'Gradient'; % [Uniform] / Gradient
    
    % --- Plotting parameters --- %
    % You can adjust this to make the plots show on your screen correctly.
    in.plotSizing = [0.1, 0.05, 0.8, 0.9];
    in.showPlots = true;
    
end

%% Call the solver!
PDEsolver

% Done.


