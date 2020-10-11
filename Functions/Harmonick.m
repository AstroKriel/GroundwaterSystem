function kCV = Harmonick(kP, ~, cfg)
    %{
        HARMONIC calculates the arithmetic average for each node face.
        
        Inputs:
            kP      relative permiability                                   (vector)
            ~
            cfg     configuration variable                                  (structure)

        Output:
            kCV     relative permiability at each face surounding a node    (structure)
    %}
    
    % Pull relevant variables from the configuration variable
    r = cfg.r;           % number of nodes in a row of our domain
    N = cfg.N;           % number of nodes in our system
    eff = cfg.efficient; % Efficient calculations

    % Initialise the relative permiability variables in (north, east) direction
    kn = zeros(N,1);
    ke = kn;
    for ii = 1:N
        if mod(ii,r) ~= 1 % north node
            kn(ii) = 2 * kP(ii) * kP(ii-1)/( kP(ii) + kP(ii-1) );
        end
        if ii <= N-r % east node
            ke(ii) = 2 * kP(ii) * kP(ii+r)/( kP(ii) + kP(ii+r) );
        end
    end

    % If we are not using efficient methodology then we need to calculate
    % ks and kw. If we are using efficient methodology then ke and kn are
    % needed only for the qe and qn calculations.
    if ~eff
        ks = [kn(2:end); 0];
        kw = [zeros(r,1); ke(1:N-r)];
        % Create a sturcture variable to return output
        kCV = struct('ks', ks, 'ke', ke, 'kn', kn, 'kw', kw);
    else
        kCV = struct('ke', ke, 'kn', kn);
    end
end
