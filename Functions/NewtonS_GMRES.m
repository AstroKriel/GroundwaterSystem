function [h, psiI, fluxI, boolConv, cfg] = NewtonS_GMRES(m, Jacobian, h, F, r, psiI, fluxI, cfg)
    % m is the number of iterations before recalculating the Jacobian using the
    % function handle given by Jacobian.
    
    % UPDATE OF NewtonS_GMRES FUNCTION
    %
    % NewtonS_GMRES function does not properly implement the GMRES
    % function. The issue is that the GMRES function always solves the 
    % Newton step within 1 iteration. The problem is most likely due to 
    % the implementation of the preconditioning matrix being incorrect.

    MaxIters = 15;
    k = 1; % initialise counter.
    % Set to be true. 
    boolConv = true;
    
    % Calculate initial function to pass to Jacobian. 
    FVM0 = F(h);
    
    % Tolerances
    tola = 1e-3;
    tolr = 1e-4;

    % Set up errors.
    err = Inf; errold = norm(FVM0,2); 
    tol=tola + tolr * errold; % Scale relative tolerance by initial error

    % Convergence rates.
    % (Testing at quadratic convergence rate)
    rho = 0.5;
    rhoT = 0;
    
    % Set hI
    hI = h;

    % % Newton iterations
    % fprintf('Newton-Shamanskii Method (m = %i)\n', m)
    % fprintf('%g: ||F(h)|| = %1.4e \n', k, errold); 

    % dummy variable
    if cfg.dt < 0.5
        linesearching = true;
        m = 1;
    else
        linesearching = false;
    end

    while k < MaxIters && err > tol
        % To Do: Set up counter other than k, so it updates every 4 not
        % mod(4).
        if mod(k, m) == 1 || rhoT > rho || k == 1
            fprintf('New Jacobian calculating at %g iterations, rhoT = %.3f\n', k, rhoT)
            J = Jacobian(F, h, FVM0, r);
            setup.type = 'crout';
            setup.milu = 'row';
            [L,U] = ilu(J,setup);
        end
 
        % Store old J, L, U, FVM0
        J_I = J; L_I = L; U_I = U; ...
        FVM0_I = FVM0;
        
        % Perform GMRES
        % Mitch: need to check all of this stuff because (according to Dale) it is ALWAYS solving after
        % one iteration which seems slightly suss.

        p = 0;           % Create counter
        dh = 0;          % Initialse Newton Step
        FVM0 = -FVM0;    % Make FVM0 negative
       
        while length(dh) == 1
            p = p + 1;

            if p == 1
                [dh] = pGMRESv02(J, FVM0, L, U); 

            elseif p == 2
                % Use old F,FVM0,L,U and h values to determine M
                [dh] = pGMRESv02(J_I,FVM0_I,L_I,U_I);
                    
            elseif p >= 3
                error('GMRES Failed to Converge')
            end
        end
        
        
        % LINE SEARCHING
        if linesearching
            lambda = 1;
            alpha = 1e-4;
            h_new = h + lambda * dh;
            while norm(F(h_new),2)^2 >= (1-2*alpha*lambda) * norm(F(h),2)^2
                sig = 0.5; % choose sig \in [0.1, 0.5]
                lambda = lambda * sig;
                h_new = h + lambda * dh;
            end
            h = h_new;
        else
            % BEFORE LINE SEARCHING
            h = h + dh;
        end

        % Calculate new F
        [FVM0, psi, flux] = F(h);
        
        % Get error. 
        err = norm(FVM0,2); 

        % Add to counter
        k = k + 1;
        
        % Makeshift adaptive timestep
        if k > 5 % only so that we don't spam the stuff.
            fprintf('%g: ||F(h)|| = %1.4e \n', k, err);
        end
        
        % Tolerance.
        rhoT = err/errold;
        errold = err;

        % Do error calculation 
        if err>tol && k == MaxIters
            boolConv = false;
        end
    end
    
    if k == 5
        cfg.dt = max(cfg.dt*0.8, cfg.dt_min); % min stepsize of a hour
        if cfg.dt == cfg.dt_min
            fprintf('\n Minimum time step reached, dt = %.2f\n', cfg.dt)
        else
            fprintf('\nAdaptive time step decreased for k = %g, dt = %.2f\n', k, cfg.dt)
        end
    elseif k == 8
        cfg.dt = max(cfg.dt*0.6, cfg.dt_min);
        fprintf('\nAdaptive time step decreased for k = %g,dt = %.2f\n', k, cfg.dt)
    end
    
    % Increase dt if quick convergence.
    if k <= 3
        if k <= 1
            cfg.dt = cfg.dt * 1.5;
        else
            cfg.dt = cfg.dt * 1.1;
        end
        
        if cfg.dt > cfg.dt_max
            cfg.dt = cfg.dt_max;
        else
            fprintf('\nAdaptive time step increased for k = %g, dt = %.2f\n', k, cfg.dt)
        end
    end
    
    % Save the convergance
    cfg.k = k;
    
    % Save the outputs
    if boolConv % we converged
        psiI = psi;
        fluxI = flux;
    else % No convergance.
        h = hI; % Return the old h. 
        % Return the old psiI too.
    end
        
end
