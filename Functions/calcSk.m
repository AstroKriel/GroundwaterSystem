function [S, kP] = calcSk(h, cfg)
    %{
        CALCSK calculates the relative saturation (S) and relative permiability (k) 
        for each node. These output vectors are predominantly dependant on h and the 
        materials around the node.

        Inouts:
            h       pressure head           (vector)
            cfg     configuration variable  (structure)

        Outputs:
            S       relative saturation     (vector [N x 4])
            kP      relative permiability   (vector)
    
    %}
    
    % Pull relevant variables from the configuration (cfg) variable
    NM      = cfg.NM;
    N       = cfg.N;
    alpha   = cfg.alpha;
    DV      = cfg.DV;
    n       = cfg.n;
    m       = cfg.m;

    % Initialise the relative saturation for all nodes
    S = zeros(size(NM));
    for j = 1:N % Loop over all nodes
        for k = 1:4 % Sub control volumes
            if ~isnan(NM(j,k)) % Ignore boundary (marked by NaN) nodes.
                if h(j) < 0 % Check head
                    S(j,k) = ( 1 + (-alpha(NM(j,k))*h(j) )^( n(NM(j,k))))^( -m(NM(j,k)) );
                else
                    S(j,k) = 1;
                end % head check
            end
        end
    end
    
    % Initialise k_all
    k_all = zeros(size(NM));

    % The k is dependent upon the 4 materials surrounding so a similar function
    % is used as per above. 
    for j = 1:N % Loop over all nodes
        for k = 1:4 % Sub control volumes
            if ~isnan(NM(j,k)) % Ignore boundary (marked by NaN) nodes.
                if h(j) < 0 % Check the head
                    k_all(j,k) = sqrt(S(j,k))*( (1-( 1-S(j,k)^( 1/m(NM(j,k)) ) )^m( NM(j,k)) )^2 );
                else
                    k_all(j,k) = 1; 
                end
            end
        end
    end

    % Calculate kP
    kP = sum(k_all.*DV,2) ./ sum(DV,2);
end


