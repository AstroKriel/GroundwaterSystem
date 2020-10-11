function kCV = Upwinding(kP, h, cfg)
    %{
        UPWINDING calculates the relative permeability for each control
        volume face, where we consider whether efficient methodology is
        requested. 

        Efficient methodology exploits the conservation of mass due to flow
        in one direction through a CV face is the negative of that flow
        from the opposite orientation.
        
        Inputs:
            kP      relative permiability                                   (vector [N x 1])
            h       pressure head                                           (vector [N x 1])
            cfg     configuration variable                                  (structure)

        Output:
            kCV     relative permiability at each face surounding a node    (structure)
    %}
    
    % Pull relevant variables from the configuration variable
    r = cfg.r; % number of nodes in a row of our domain
    N = cfg.N; % number of nodes in our system
    Z = cfg.Z; % Z values
    eff = cfg.efficient; % efficient calculations?
    
    H = vec2mat(h,r)' + Z;
    
    % Now we want to take the H values that are bigger than the east
    % neighbours
    
    % Initialise the relative permiability variables in (north, east) direction
    kn = zeros(N,1);
    ke = kn;
    
    for ii = 1:N 
        
        if mod(ii,r) ~= 1 % if == 1, then on the north boundary
            if H(ii-1) > H(ii)
                kn(ii) = kP(ii-1);
            else
                kn(ii) = kP(ii);
            end
        end
        
        if ii <= N-r % if > N-r on the east boundary.
            if H(ii+r) > H(ii)
                ke(ii) = kP(ii+r);
            else
                ke(ii) = kP(ii);
            end
        end
        
    end

    % If we are using the efficient code methodology we calculate the ks
    % and kw control volume face values as well. We then save these values
    % in a struct to return.
    
    if ~eff
        ks = [kn(2:end); 0];
        kw = [zeros(r,1); ke(1:N-r)];
        kCV = struct('ks', ks, 'ke', ke, 'kn', kn, 'kw', kw);
    else
        kCV = struct('ke', ke, 'kn', kn);
    end
    
end
