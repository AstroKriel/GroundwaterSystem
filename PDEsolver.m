% Are we loading from a file?
if load_OLD
    % Get ready to start computations again.
    % Mitch: check the logic flow for this

    % Let QI = Q;
    % Let psi = psiI
    QI = Q; % Assume momentary steady state?
    psi = psiI;
else
    % Produces a valid configuration based on the inputs set in the user's Master.
    % Contains validation logic to ensure that things are done correctly.
    [fm, cfg] = checkCFG(in, cfgMesh);
    clear in
    
    % Initial calculations
    [cfg, psi, S, kP, h, Q, fluxI] = InitialCalculations(cfg);

    % Unpack and Initialise
    psiI = psi;
    QI = Q;
    t_total = 0;
    cfg.t_total = t_total;
end

% Unpack some variables
t_max = cfg.t_max;
r = cfg.r;
N = cfg.N;

% Initialise vectors that track the performance/state of the system
% Mitch: this will tie into analysis
psi_vec      = NaN(min(t_max,1e4), 1); % track the average water content in the system
psi_ini      = sum(psi .* cfg.CV) / (cfg.L1 * cfg.L2);
psi_vec(1)   = psi_ini;
t_vec        = NaN(min(t_max,1e4), 1); % track the elapsed time
t_vec(1)     = 0;
dt_vec       = NaN(min(t_max,1e4), 1); % track the change in dt
dt_vec(1)    = cfg.dt;
k_vec        = NaN(min(t_max,1e4), 1); % track the change in steps req for convergence in k
k_vec(1)     = 0;
loopTimings  = NaN(min(t_max,1e4), 1); % track timing
river_vec    = NaN(min(t_max,1e4), 1); % track average river flux
river_vec(1) = cfg.dt * mean( fluxI(1:cfg.Nriver, 2)  .* cfg.CV(1:cfg.Nriver) );
rain_vec     = NaN(min(t_max,1e4), 1); % track average rain flux
rain_vec(1)  = cfg.dt * mean( fluxI(1:r:end, 1) .* cfg.CV(1:r:end) );
csg_vec      = NaN(min(t_max,1e4), 1); % track average csg flux
csg_vec(1)   = cfg.dt * mean( fluxI(cfg.NCSG:end, 4) .* cfg.CV(cfg.NCSG:end) );
Q_vec        = NaN(min(t_max,1e4),1); % track Q flow
Q_vec(1)     = cfg.dt * mean( Q(1:r:end) .* cfg.CV(1:r:end) );
tot_vec      = NaN(min(t_max,1e4),1); % track total flow


%% Are we displaying plots?
if cfg.showPlots
    set(gcf, 'Units', 'Normalized', 'OuterPosition', cfg.plotSizing);

    % figure scales depending on the screensize and dimensions. The figure starts (bottom left)
    % at 10% the width of screen and 20% the height of screen, with the dimensions of the screen 
    % defined to be 70% the width and 60% the height of the screen.

    %% Initialise subplots
    % Mitch: add in cfg variable allowing the plots to be toggled on/off for speed reasons. Done.
    h1 = subplot(3,2,1);
    plot(t_vec, psi_vec, '-', 'LineWidth', 1.5)
    title('Average Water Content', ...
            'Interpreter', 'LateX', 'FontSize', 20)
        hold on
        psi_cur = psi_ini + cfg.t_total * cfg.rainfallSurface / cfg.L2;
        plot([0, cfg.t_total], [psi_ini, psi_cur], '--', 'LineWidth', 1.5)
        hold off
        legend('System''s Average Water Content', 'Predicted Water Content in Simple System', ...
            'Location', 'SouthEast')    
        xlabel('Time (days)', 'Interpreter', 'LateX', 'FontSize', 16)
        ylabel('Average Water Content ($m^2$)', 'Interpreter', 'LateX', 'FontSize', 16)
        
    h2 = subplot(3,2,2);
    title('Change in Step Size and k', 'Interpreter', 'LateX', 'FontSize', 20)
        xlabel('Iterations', 'Interpreter', 'LateX', 'FontSize', 16)
        % left axis
        yyaxis left
        plot(dt_vec, '-', 'LineWidth', 1.5)
        ylabel('Step Size (dt)', 'Interpreter', 'LateX', 'FontSize', 16)
        yyaxis right
        plot(1:length(k_vec), k_vec, '--', 'LineWidth', 1.5)
        ylim([0, 15])
        ylabel('Convergence of k', 'Interpreter', 'LateX', 'FontSize', 16)

    h3 = subplot(3,2,3);
    surf(cfg.X, cfg.Z, vec2mat(h,r)', 'EdgeColor', 'none', 'FaceColor', 'interp')
    view(2)
    title('Pressure Head', 'Interpreter', 'LateX', 'FontSize', 20)
        xlabel('Width (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
        ylabel('Height (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
        colorbar
        caxis([-15 0])

    h4 = subplot(3,2,4);
    surf(cfg.X, cfg.Z, vec2mat(psiI,r)', 'EdgeColor', 'none', 'FaceColor', 'interp')
    view(2)
    title('Water Content', 'Interpreter', 'LateX', 'FontSize', 20)
        xlabel('Width (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
        ylabel('Height (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
        colorbar
        caxis([0 0.4686])

    h5 = subplot(3,2,5);
    plot(t_vec, river_vec, '-', 'LineWidth', 1.5)
        hold on
        plot(t_vec, rain_vec, '--', 'LineWidth', 1.5)
    title('Water Entering/Leaving System Across Boundaries', ...
            'Interpreter', 'LateX', 'FontSize', 20)
        xlabel('Time (days)', 'Interpreter', 'LateX', 'FontSize', 16)
        ylabel('Water Content', 'Interpreter', 'LateX', 'FontSize', 16)
        legend('River Flux', 'Rain Flux', ...
            'Location', 'SouthWest')

    h6 = subplot(3,2,6);
    surf(cfg.X, cfg.Z, vec2mat(psiI./cfg.PsatNode,r)', 'EdgeColor', 'none', 'FaceColor', 'interp')
    view(2)
    % surf(cfg.X,cfg.Z,vec2mat((psiI./cfg.PsatNode),r)')
    title('$\psi/\psi_{\mbox{sat}}$', 'Interpreter', 'LateX', 'FontSize', 20)
        xlabel('Width (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
        ylabel('Height (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
        colorbar
        caxis([0 1])

end

%% Simulation
iter = cfg.iter;

while cfg.t_total < t_max % Time Loop (converted from years to days)
    
    % update iteration varibale
    iter = iter + 1;

    % Flux (turned on and off by csgFluxTime etc.)
    if cfg.t_total > cfg.csgFluxTime(2) || cfg.t_total < cfg.csgFluxTime(1)% turn off
        cfg.csgFlux = false;
    elseif cfg.t_total > cfg.csgFluxTime(1) && ... turn on
                all(psiI(cfg.NCSG:N) ./ cfg.PsatNode(cfg.NCSG:N) > 0.2)
        cfg.csgFlux = true;
    end
    
    % Pump
    if cfg.t_total > cfg.pumpOnTime(2) % turn off
        cfg.pumpOn = false;
    elseif cfg.t_total > cfg.pumpOnTime(1) % add saturation condition?
        cfg.pumpOn = true;
    end
    
    % Rain
    if cfg.t_total > cfg.rainFluxTime(2) % turn off
        cfg.rainFlux = false;
    elseif cfg.t_total > cfg.rainFluxTime(1) % turn on
        cfg.rainFlux = true;
    end
        
    % --- SOLVE --- %
    tic
    
    % --- calculate daily rainfall --- %
    cfg.rainfallSurface = cfg.rainFunc(cfg.t_total, cfg);
    
    % Track total ellapsed time
    cfg.t_total = cfg.t_total + cfg.dt;
    
    % Because fluxI, QI and psiI change each loop, we update F(h). 
    F = @(h) FVM(h, Q, QI, fluxI, psiI, cfg);
    
    % Newton solve: pass Jacobian, and new F handle. 
    % Set boolConv to be false to initialise
    boolConv = false;
    tic
    while ~boolConv
        
        [h, psiI, fluxI, boolConv, cfg] = NewtonS(4, cfg.Jacobian, h, F, r, ...
                                                      psiI, fluxI, cfg);
        
        if ~boolConv
            cfg.dt = cfg.dt/2; % Big cut 
            if cfg.dt < cfg.dt_min
                error('Did not converge')
            end
        end
    end
    ct = toc;
    
    % Q is called every loop is evapo is enabled.
    if cfg.evapo
        QI = Q; % Save old Q
        Q = Qcalc(psiI, cfg); % Calculate new Q
    end
    
    % --- ANALYSIS --- %
    % Mitch: also need to toggle this code with a config variable.
    % Mitch: ideally, should it save the analysis results etc with the save points? Check.

    % Look at the amount of water coming in/leaving our system
    river_vec(iter + 1) = cfg.dt * mean( fluxI(1:cfg.Nriver,4)  .* cfg.CV(1:cfg.Nriver) );
    rain_vec(iter + 1)  = - cfg.dt * mean( fluxI(1:r:end,3) .* cfg.CV(1:r:end) );
    csg_vec(iter + 1)   = cfg.dt * mean( fluxI(cfg.NCSG:end,4) .* cfg.CV(cfg.NCSG:end) );
    Q_vec(iter + 1)   = cfg.dt * mean( Q(1:r:end) .* cfg.CV(1:r:end) );
    tot_vec = Q_vec + rain_vec + river_vec + csg_vec;
    
    % Looking at the amount of average water content in our system
    psi_vec(iter + 1)   = sum(psiI.*cfg.CV) / (cfg.L1*cfg.L2);
    
    % Tracking other variables
    t_vec(iter + 1)     = cfg.t_total; % time
    dt_vec(iter + 1)    = cfg.dt;      % change in step size
    k_vec(iter + 1)     = cfg.k;       % reltive permeability

    % Commenting out starts here
    if cfg.showPlots
        % --- PLOTS --- %
        % plot the average water content (expected vs calculated)
        subplot(3, 2, 1)
        hold on
        cla
        plot(t_vec, psi_vec, '-', 'LineWidth', 1.5)
        psi_cur = psi_ini + cfg.t_total * cfg.rainfallSurface / cfg.L2;
        plot([0, cfg.t_total], [psi_ini, psi_cur], '--', 'LineWidth', 1.5)
        hold off
        
        % plot the change in step size
        subplot(3, 2, 2)
        % title('Change in Step Size and k', 'Interpreter', 'LateX', 'FontSize', 20)
        % xlabel('Iterations', 'Interpreter', 'LateX', 'FontSize', 16)
        % left axis
        yyaxis left
        hold on
        cla
        plot(dt_vec, '-', 'LineWidth', 1.5)
        ylim([0, 1.1*max(dt_vec)])
        % ylabel('Step Size (dt)', 'Interpreter', 'LateX', 'FontSize', 16)
        % right axis
        yyaxis right
        hold on
        cla
        plot(1:length(k_vec), k_vec, '--', 'LineWidth', 1.5)
        % ylim([0, 15])
        % ylabel('Convergence of k', 'Interpreter', 'LateX', 'FontSize', 16)
        
        % plot the pressure head in the system
        subplot(3, 2, 3, h3)
        hold on
        cla
        surf(cfg.X, cfg.Z, vec2mat(h,r)', 'EdgeColor', 'none') % , 'FaceColor', 'interp'
        view(2)
        % title('Pressure Head', 'Interpreter', 'LateX', 'FontSize', 20)
        % xlabel('Width (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
        % ylabel('Height (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
        
        % plot the water content in the system
        subplot(3, 2, 4, h4)
        hold on
        cla
        surf(cfg.X, cfg.Z, vec2mat(psiI,r)', 'EdgeColor', 'none', 'FaceColor', 'interp')
        view(2)
        % title('Water Content', 'Interpreter', 'LateX', 'FontSize', 20)
        % xlabel('Width (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
        % ylabel('Height (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
        % colorbar
        
        subplot(3, 2, 5, h5)
        hold on
        cla
        plot(t_vec, cfg.dt .* river_vec, '-', 'LineWidth', 1.5)
        hold on
        plot(t_vec, cfg.dt .* rain_vec, '--', 'LineWidth', 1.5)
        % title('Water Entering/Leaving System Across Boundaries', ...
        %       'Interpreter', 'LateX', 'FontSize', 20)
        % xlabel('Time (days)', 'Interpreter', 'LateX', 'FontSize', 16)
        % ylabel('Water Content', 'Interpreter', 'LateX', 'FontSize', 16)
        % legend('River Flux', 'Rain Flux', ...
        %        'Location', 'SouthWest')
        hold off
        
        % plot the water content in the system
        subplot(3, 2, 6, h6)
        hold on
        cla
        surf(cfg.X, cfg.Z, vec2mat(psiI ./ cfg.PsatNode,r)', 'EdgeColor', 'none', 'FaceColor', 'interp')
        view(2)
        % surf(cfg.X,cfg.Z,vec2mat((psiI./cfg.PsatNode),r)')
        % title('$\psi/\psi_{\mbox{sat}}$', 'Interpreter', 'LateX', 'FontSize', 20)
        % xlabel('Width (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
        % ylabel('Height (meters)', 'Interpreter', 'LaTeX', 'FontSize', 16)
        % colorbar
        
        drawnow
    end
    
    loopTimings(iter) = toc;
    % Average Time (AT), Current Time (CT)
    % fprintf('AT = %.4f \t CT = %.4f\n',...
    %         mean(loopTimings, 'omitnan'), loopTimings(iter)) 
    
    % SAVE
    if mod(iter,cfg.saveITER) == 0
        cfg.iter = iter;
        cfg = SaveData(h, psiI, Q, fluxI, cfg, fm);
    end
    
    if mod(iter, 10) == 0 && cfg.saveH
        save_h(cfg, h, cfg.t_total)
    end
    
    % Mitch: can we add any extra useful output without spamming?
    fprintf('Current time: %g\t\tLoopTime: %g\n', cfg.t_total, ct);
    
end


% Save data
cfg.iter = iter;
cfg = SaveData(h, psiI, Q, fluxI, cfg, fm);