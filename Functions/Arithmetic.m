function kCV = Arithmetic(kP, ~, cfg)
    %{
        ARITHMATIC calculates the arithmetic average for each node face.
        
        Inputs:
            kP      relative permeability                                   (vector)
            ~
            cfg     configuration variable                                  (structure)

        Output:
            kCV     relative permeability at each face surrounding a node   (structure)
    %}
    
    % Pull relevant variables from the configuration variable
    r = cfg.r; % number of nodes in a row of our domain
    N = cfg.N; % number of nodes in our system

    % Initialise the relative permiability variables in (north, east) direction
    kn = zeros(N,1);
    ke = kn;
    for ii = 1:N
        if mod(ii,r) ~= 1 % north node
            kn(ii) = (kP(ii) + kP(ii-1))/2;
        end
        if ii <= N-r % east node
            ke(ii) = (kP(ii)+ kP(ii+r))/2;
        end
    end

    % create the south and west nodes
    ks = [kn(2:end); 0];
    kw = [zeros(r,1); ke(1:N-r)];

    % create a sturcture variable to return output
    kCV = struct('ks', ks, 'ke', ke, 'kn', kn, 'kw', kw);
end
