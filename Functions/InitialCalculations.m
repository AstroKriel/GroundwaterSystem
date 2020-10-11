function [cfg, psi, S, kP, h, Q, fluxI] = InitialCalculations(cfg)

    Z = cfg.Z;
    N = cfg.N;
    rainFunc = cfg.rainFunc;

    % --- Initial Conditions --- %

    % pressure head at top and bottom. 
    h_bot = -5; h_top = -10;

    % Calculate the initial conditions for all nodes.
    h = initial_h_conditions(h_bot, h_top, Z, cfg.L2);

    % Clear unecessary variables from workspace.
    clear h_bot h_top

    % Initialise S, kp
    [S, kP] = calcSk(h, cfg);

    % More initialisation stuff
    kCV = cfg.Solvek(kP, h, cfg);

    % Initialise the Water Content (PSI) at each node 
    psi = calcPsi(h, S, cfg);

    % Calculate the Psat of each node relative to surrounding materials
    cfg.PsatNode = zeros(size(psi));
    for iter = 1:N % Loop over all nodes
        counter = 0;
        for k = 1:4 % All materials
            if ~isnan(cfg.NM(iter,k))
                counter = counter + 1;  % Add one to counter for average
                cfg.PsatNode(iter) = cfg.PsatNode(iter) + cfg.Psat(cfg.NM(iter,k)); % Add saturation
            end
        end
        cfg.PsatNode(iter) = cfg.PsatNode(iter) / counter; % Divide by counter for average
    end
    clear counter

    % Calculate rainfall at surface.
    cfg.rainfallSurface = rainFunc(0, cfg);

    % Calculate yearly rain.
    cfg.yearRain = cfg.rainfallSurface;

    % This is correct because the line above
    for t = 1:364
        cfg.yearRain = cfg.yearRain + rainFunc(t, cfg);
    end

    % Calculate the initial conditions for Flow Rate.
    if cfg.evapo
        Q = Qcalc(psi, cfg);
    else
        Q = zeros(N,1);
    end

    % Old Flux for initial calc.
    fluxI = FluxCalculator(h, kCV, 0, cfg);

end