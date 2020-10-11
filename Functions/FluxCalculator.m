function flux = FluxCalculator(h, kCV, psi, cfg)
    %{
        FLUXCALCULATOR
        Inputs:
            h       pressure head               (vector)
            kCV     
            dt      time step                   
            cfg     configuration varibale      (structure)

        Output:
            flux
    %}

    % Pull relevant variables from the configuration (cfg) variable
    N         = cfg.N;
    r         = cfg.r;
    c         = cfg.c;
    Z         = cfg.Z;
    NM0       = cfg.NM0;
    x_CSG     = cfg.x_CSG;
    H_CSG     = cfg.H_CSG;
    riverFlux = cfg.riverFlux;
    csgFlux   = cfg.csgFlux;
    Nriver    = cfg.Nriver;
    NCSG      = cfg.NCSG;
    K_R       = cfg.K_R;
    H_R       = cfg.H_R;
    R_b       = cfg.R_b;
    x_R       = cfg.x_R;
    Kxx       = cfg.Kxx;
    Kzz       = cfg.Kzz;
    delta_xw  = cfg.delta_xw;
    delta_xe  = cfg.delta_xe;
    delta_zs  = cfg.delta_zs;
    delta_zn  = cfg.delta_zn;
    efficient = cfg.efficient;
    rainflux  = cfg.rainFlux;
    pumpOn    = cfg.pumpOn;
    f_R       = cfg.f_R; 
    pumpNu    = cfg.pumpNu;
    pumpNl    = cfg.pumpNl;

    % Get hydraulic head. 
    H = h(:) + Z(:);

    % NORTH
    % Rain flux is negative (i.e. opposite direction to positive orientation).
    qn = zeros(size(Z)); % initialise
    if rainflux
        % The nodes at z = L2 are all boundary nodes where the north flux is given by the rainfall.
        drain = cfg.rainfallSurface; % This is total rainfall over all nodes (note independent of time).
        xCV   = delta_xe(1,:) + delta_xw(1,:); % get CV width in x

        % Check whether or not the node is saturated.
        % If it is we set the rainfall to 0 for that node.
        % Else nothing happens.
        psi = vec2mat(psi, r)';
        PsatNode = vec2mat(cfg.PsatNode, r)';
        boolSat = psi(1,:) >= PsatNode(1,:);
        if any(boolSat)
            xCV(boolSat) = 0;
        end
        qn(1,:) = -( xCV * drain ); % Scale by the control volume width.
    end

    % EAST
    % Initialise matrix;
    qe = zeros(N,1); 

    % CSG flux. 
    if csgFlux
        q_csg = Kxx(NM0(N,1)) * ( H_CSG - H(NCSG:end) ) / x_CSG;
        q_csg(1) = q_csg(1) * ( delta_zn(NCSG) + delta_zs(NCSG) ) / delta_zs(NCSG);
        qe(NCSG:N) = q_csg;
    end

    % West and South.
    if ~efficient
        qw = qe;
        qs = qe;
    end

    % Now we loop over all nodes and calculate the rest of the flux.
    for j = 1:N
        % NORTH
        if mod(j,r) ~= 1
            qn(j) = -kCV.kn(j) * ...
                ( ...	% Open bracket
                Kzz(NM0(j,1)) * delta_xw(j) + Kzz(NM0(j,2)) * delta_xe(j) ...
                 ) ...  % Close bracket
                 * ( H(j-1)-H(j) )/( delta_zn(j) + delta_zs(j-1) );
        end

        % SOUTH
        % Run this if efficiency is off.
        if mod(j,r) ~= 0 && ~efficient 
            qs(j) = -kCV.ks(j) * ...
                (...   % Open bracket
                Kzz(NM0(j,3)) * delta_xw(j) + Kzz(NM0(j,4)) * delta_xe(j) ...
                ) ...  % Close bracket
                * ( H(j)-H(j+1) )/( delta_zs(j) + delta_zn(j+1) );
        end

        % EAST
        if j <= N-r 
            % internal nodes. 
            qe(j) = -kCV.ke(j) * ...
                ( ...  % Open bracket
                Kxx(NM0(j,2)) * delta_zn(j) + Kxx(NM0(j,4)) * delta_zs(j) ...
                ) ...  % Close bracket
                * ( H(j+r)-H(j) )/( delta_xe(j) + delta_xw(j+r) );
        end

        % WEST
        % Run this if efficiency is off.
        if j > r && ~efficient
            qw(j) = -kCV.kw(j) * ( ...  % Open bracket
                        Kxx(NM0(j,1)) * delta_zn(j) + Kxx(NM0(j,3)) * delta_zs(j) ...
                                ) ...   % Close bracket
                      * ( H(j)-H(j-r) )/( delta_xw(j) + delta_xe(j-r) );
        elseif j <= Nriver && riverFlux && ~efficient
            % If in river region, & river flux is active & not efficient.
            % Calculate river nodes.
            qw(j) = K_R * ( H_R - max(R_b, H(j)) ) / x_R;

            % Scale by length of CV face.
            if j == 1
                %qw(j) = qw(j) / 2;
            elseif j == Nriver
                qw(j) = qw(j) * delta_zn(j) / ( delta_zn(j) + delta_zs(j) );
            end
        end
    end
    
    % Consider whether the pump is on. If so correct pumping rates. 
    

    % Calculate qs
    if efficient
        qs = [qn(2:end,:); zeros(1,c)];
    end
    
    if pumpOn
        qn(pumpNl:pumpNu)       = f_R .* qn(pumpNl:pumpNu);
        qn(pumpNl-r:pumpNu-r)   = f_R .* qn(pumpNl-r:pumpNu-r);
        qn(pumpNl+r:pumpNu+r)   = f_R .* qn(pumpNl+r:pumpNu+r);
        qs(pumpNl:pumpNu)       = f_R .* qs(pumpNl:pumpNu);
        qs(pumpNl-r:pumpNu-r)   = f_R .* qs(pumpNl-r:pumpNu-r);
        qs(pumpNl+r:pumpNu+r)   = f_R .* qs(pumpNl+r:pumpNu+r);
    end

    if riverFlux && efficient
        % Calculate qw
        q_river = K_R * (H_R - max(ones(Nriver,1) .* R_b, H(1:Nriver)))/x_R;

        % Scale the ones on the edge between no flux and flux.
        q_river(end) = q_river(end) * delta_zn(Nriver)/ ( delta_zn(Nriver) + delta_zs(Nriver) );

        qw = [q_river; zeros(r - Nriver,1); qe(1:N-r)];
    elseif efficient
        qw = [zeros(r,1); qe(1:N-r)];
    end

    % Thus we now have the flux for each node. Vectorise North and South.
    flux = [qs(:), qe, qn(:), qw];
     
end
