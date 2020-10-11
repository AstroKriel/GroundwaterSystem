function psi = calcPsi(h, S, cfg)
    %{
        CALCPSI calculates all the water content for all the nodes in our system
    %}

    % Pull relevant variables from the configuration (cfg) variable
    N = cfg.N;
    DV = cfg.DV;
    Psat = cfg.Psat;
    Pres = cfg.Pres;
    NM = cfg.NM;

    % Initialise the water content vector
    psi_all = zeros(size(NM));
    for j = 1:N % Loop over all nodes
        for k = 1:4 % Sub control volumes
            if ~isnan(NM(j,k)) % Ignore boundary (marked by NaN) nodes.
                if h(j) < 0 % Check head
                    psi_all(j,k) = Pres(NM(j,k)) + S(j,k) * (Psat(NM(j,k)) - Pres(NM(j,k)));
                else
                    psi_all(j,k) = Psat(NM(j,k));
                end
            end
        end
    end

    % Scale (?)
    psi = sum(psi_all.*DV,2) ./ sum(DV,2);
end
