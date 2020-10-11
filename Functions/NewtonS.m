function [h, psiI, fluxI, boolConv, cfg] = NewtonS(m, Jacobian, h, F, r, psiI, fluxI, cfg)
    % m is the number of iterations before recalculating the Jacobian using the
    % function handle given by Jacobian.

    k = 1;           % initialise counter.
    boolConv = true; % Set to be true.
    MaxIters = 15;    
    
    % Calculate initial function to pass to Jacobian. 
    FVM0 = F(h);
    
    % Tolerances
    tola = 1e-1;
    tolr = 1e-2;

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

    % Mitch: this makes sense; assuming it's something we can write about? or
    % is it too trivial to worry about / assumed that we'll be doing this.
    % Only linesearching if dt < 0.5
    if cfg.dt < 0.5
        linesearching = true;
        m = 1;
    else
        linesearching = false;
    end

    while k < MaxIters && err > tol

        if mod(k, m) == 1 || rhoT > rho || k == 1
            % Optional debugging
            % fprintf('New Jacobian calculating at %g iterations, rhoT = %.3f\n', k, rhoT)
            J = Jacobian(F, h, FVM0, r);
        end
        
        % Apply the Jacobian
        dh = J\(-FVM0);
        
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
        % k > 5 prevents cluttering console output; only shows output
        % in exceptional circumstances (ie, more than 5 loops)
        if k > 5
            fprintf('%g: ||F(h)|| = %1.4e \n', k, err);
        end
        
        % Tolerance.
        rhoT = err / errold;
        errold = err;

        % Do error calculation 
        if err > tol && k == MaxIters
            boolConv = false;
        end
    end
    
    if k == 5
        cfg.dt = max(cfg.dt*0.8, cfg.dt_min); % min stepsize of a hour
        if cfg.dt == cfg.dt_min
            fprintf('\nMinimum time step reached, dt = %.2f\n', cfg.dt)
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
        if exist('psi',1)
            error('Calculations not correct, check h values')
        end
        
        psiI = psi;
        fluxI = flux;
    else % No convergance.
        h = hI; % Return the old h. 
        % Return the old psiI too.
    end
        
end
