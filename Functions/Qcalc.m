function Q = Qcalc(psi, cfg)
    %{
        Qcalc calculates flow ...

        Evapotranspiration.
        The regions in which evapotranspiration occur:
        repairing vegetation:   x \in [50, 100]
                                z \in [85, 100]
        crop zone:              x \in [100, 300]
                                z \in [95, 100]
        pine plantation:        x \in [300, 500]
                                z \in [90, 100]

        Inouts:
            psi     water content           (vector)
            cfg     configuration variable  (structure)

        Output:
            Q       Flow                    (vector)

    %}

    % Q calc is only called once so optimization isn't worth much at all here.
    % Mitch: Need to double-check that the above is true - checked, was false.
    % This function is called periodically as we move forward in time.

    % Pull relevant variables from the configuration (cfg) variable
    x               = cfg.x;
    z               = cfg.z;
    ell1            = cfg.ell1;
    ell2            = cfg.ell2;
    ell3            = cfg.ell3;
    L2              = cfg.L2;
    R_1p            = cfg.R_1p;
    R_2p            = cfg.R_2p;
    R_3p            = cfg.R_3p;
    pumpOn          = cfg.pumpOn;
    yearRain        = cfg.yearRain; % Created each year from rainFunc
    rainfallSurface = cfg.rainfallSurface;
    PsatNode        = cfg.PsatNode;
    N               = cfg.N;
    P_R             = cfg.P_R;
    pumpNu          = cfg.pumpNu;
    pumpNl          = cfg.pumpNl;
    CV              = cfg.CV;

    % The percentages.
    R_1 = R_1p * rainfallSurface;
    R_2 = R_2p * rainfallSurface;
    R_3 = R_3p * rainfallSurface;

    % Size of the Flow Matrix.
    Q = zeros(length(z), length(x));

    % Loop over every node in the system
    % Mitch: we can optimize this by avoiding looping the lower half 
    % of the system which is not affected.
    for jx = 1:1:cfg.c
        for jz = 1:1:cfg.r
            % Does the control volume around this node collide with any
            % boundaries? Get the relative amount of each zone type in
            % the volume for this node.
            [vegPV, cropPV, pinePV] = getPercentageVolAroundNode(jx, jz, cfg);
            
            vegQ  = - R_1 * ( z(jz) - L2 + ell1 )^2 / ell1^2;
            cropQ = - R_2 * ( z(jz) - L2 + ell2 )^2 / ell2^2;
            pineQ = - R_3 * ( z(jz) - L2 + ell3 )^2 / ell3^2;
            
            % Add relative amounts of Q's together to get Q at this node.
            % If there is no veg, crop, or pine, this will evaluate to zero.
            Q(jz, jx) = vegPV * vegQ + ...
                        cropPV * cropQ + ...
                        pinePV * pineQ;
            
        end
    end

    % See Page 6 of the Project Description
    for jj = 1:N
        if psi(jj)/PsatNode(jj) <= 0.5
            Q(jj) = 0;
        end
    end

    if pumpOn
        % dt is measured in days. thus we take a total amount of water given by
        % the following:
        
        % This calculates the percentage of rainfall we take based off the
        % yearly values of the rainfall. 
        pumpTotal = - P_R * yearRain / 365; % Convert to correct units
        fprintf('pumpTotal = %g\n', pumpTotal); 

        % This is spread evenly over the nodes we pump from.
        % Find the coordinates.
        Q(pumpNl:pumpNu) = CV(pumpNl:pumpNu) * pumpTotal / abs((pumpNu - pumpNl));
    end

    % Vectorise output
    Q = Q(:);
end

function [veg, crop, pine] = getPercentageVolAroundNode(jx, jz, cfg)
    % Author: Mitchell Johnson
    % Computes the relative amount (between 0 and 1) that each type of region
    % occupies within the control volume at jx, jz in our system.
    % Eg, veg = 1, crop = 0, pine = 0 means that the region is entirely within
    % the vegetation area.

    % Mini Diagram:
    % --------------------------------------------------
    %     |   1   |------------------------|     3     |
    %     |       |           2            |-----------|
    %     |-------|

    % Unpack Variables
    x               = cfg.x;
    z               = cfg.z;
    ell1            = cfg.ell1;
    ell2            = cfg.ell2;
    ell3            = cfg.ell3;
    L2              = cfg.L2;
    r               = cfg.r;

    % Get the volume and delta sizes around this node
    dxe = cfg.delta_xe(jz, jx);
    dxw = cfg.delta_xw(jz, jx);
    dzn = cfg.delta_zn(jz, jx);
    dzs = cfg.delta_zs(jz, jx);
    vol = cfg.CV((jx - 1) * r + jz); % pre-computed

    % THIS IS CORRECT I HAVE DOUBLE CHECKED
    % Compute veg
    vegOverlap = overlapArea(x(jx) - dxw, z(jz) - dzs, x(jx) + dxe, z(jz) + dzn, ...
                             50, L2 - ell1, 100, L2); % veg region
    veg = vegOverlap / vol;

    % Compute crop
    cropOverlap = overlapArea(x(jx) - dxw, z(jz) - dzs, x(jx) + dxe, z(jz) + dzn, ...
                              100, L2 - ell2, 300, L2); % crop region
    crop = cropOverlap / vol;

    % Compute pine
    pineOverlap = overlapArea(x(jx) - dxw, z(jz) - dzs, x(jx) + dxe, z(jz) + dzn, ...
                              300, L2 - ell3, 500, L2); % pine region
    pine = pineOverlap / vol;
end

function area = overlapArea(x1, y1, x2, y2, x1s, y1s, x2s, y2s)
    % Author: Mitchell Johnson
    % Computes the overlap between two rectangles given by x1, y1, x2, y2, and 
    % x1s, y1s, x2s, y2s respectively.

    x_overlap = max(0, min(x2, x2s) - max(x1, x1s));
    y_overlap = max(0, min(y2, y2s) - max(y1, y1s));
    area = x_overlap * y_overlap;
end