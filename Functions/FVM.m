function varargout = FVM(h, Q, QI, fluxI, psiI, cfg)
    %{ 
        Finite Volume Model (FVM)

        Input:
            h           head            (current)   (vector)
            Q           flow            (current)   (vector)
            QI          flow            (past)      (vector)
            fluxI       flux            (past)      (vector)
            psiI        water content   (past)      (vector)
            cfg         configuration variable      (structure)

        Output: 
            varargout - the number of outputs are dynamically determined depending on
            which configuration of outputs are requested.
            i.e.    For Jacobian:   return F only. 
                    For NewtonS:    return F, psi, and flux. 
    %}
    
    % Pull relevant variables from the configuration (cfg) variable
    Solvek   = cfg.Solvek;
    theta_1  = cfg.theta_1;
    theta_2  = cfg.theta_2;
    dt = cfg.dt;
    CV = cfg.CV;

    % Calculate S and kP for each node
    [S, kP] = calcSk(h, cfg);

    % Calculate k for each node's control volumes
    kCV = Solvek(kP, h, cfg); % Solvek caluclated using: Arithmatic.m

    % Calculate the volumetric water content at each node 
    psi = calcPsi(h, S, cfg);

    % Calculate the flux
    flux = FluxCalculator(h, kCV, psi, cfg);
    
    % flux = [qs, qe, qn, qw];
    qs = flux(:,1);
    qe = flux(:,2);
    qn = flux(:,3);
    qw = flux(:,4);

    % Old flux
    qsI = fluxI(:,1);
    qeI = fluxI(:,2);
    qnI = fluxI(:,3);
    qwI = fluxI(:,4);

    % FVM from discretisation
    F = psi - psiI + dt * ( ... % open dt bracket
                    theta_1.*( qn - qs + qe - qw )./CV - theta_2*Q ...
               + (1-theta_1).*( qnI - qsI + qeI - qwI )./CV - (1-theta_2)*QI ...
                          );    % Close dt bracket

    % For the Jacobian we only want to return F. 
    varargout{1} = F;

    % For the Newtons we return psi and flux as well. 
    if nargout == 3
        varargout{2} = psi;
        varargout{3} = flux;
    end
end
