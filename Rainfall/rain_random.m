function r_am = rain_random(t, cfg)
% This calculates the rainfall on a given day using a combination of cosine
% terms and an added random amount to simulate large amounts of rainfall 
% r_am refers to rainfall amount
% t refers to the day
% cfg is used for consistency with the other functions

% Calculate the average daily rainfall 
% This uses the least squares fitted model scaled down to account for
% downpour events
if isfield(cfg, 'droughtStart')
    if t >= cfg.droughtStart && t <= cfg.droughtStart + cfg.droughtLength % Enforce drought conditions
        if rand < 0.01 % Allow a small chance of reduced rain
            r_am = 0.5*(1.883*10^(-3) + 2.702*10^(-5)*cos(2*pi*t/1825) + 1.023*10^(-3)...
                *cos(2*pi*t/365) - 1.719*10^(-4)*cos(2*pi*t/90) - 6.189*10^(-5)*cos(2*pi*t/30));
            
            % Check to see if a downpour occurs
            if mod(t,365) <= 60 || mod(t,365) >= 335 % Summer rains
                if rand(1) > 0.91
                    r_am =  0.02513 + (-1 + (1+1)*rand(1))*0.01705847; % The downpour event occurs
                end
            elseif mod(t,365) >= 61 || mod(t,365) <= 151 % Autumn rains
                if rand(1) > 0.95
                    r_am = 0.02506 + (-1 + (1+1)*rand(1))*0.01560596; % The downpour event occurs
                end
            elseif mod(t,365) >= 152 || mod(t,365) <= 262 % Winter rains
                if rand(1) > 0.97
                    r_am = 0.02016 + (-1 + (1+1)*rand(1))*0.01085544; % The downpour event occurs
                end
            else % Spring rains
                if rand(1) > 0.93
                    r_am = 0.02162 + (-1 + (1+1)*rand(1))*0.01228265; % The downpour event occurs
                end
            end
            r_am = r_am*0.5;
        else
            r_am = 0;
        end
    else % Normal rainfall occurs
    r_am = 0.5*(1.883*10^(-3) + 2.702*10^(-5)*cos(2*pi*t/1825) + 1.023*10^(-3)...
        *cos(2*pi*t/365) - 1.719*10^(-4)*cos(2*pi*t/90) - 6.189*10^(-5)*cos(2*pi*t/30));
    
    end
else        
    % Check to see if a downpour occurs
    if mod(t,365) <= 60 || mod(t,365) >= 335 % Summer rains
        if rand(1) > 0.91
            r_am =  0.02513 + (-1 + (1+1)*rand(1))*0.01705847; % The downpour event occurs
        end
    elseif mod(t,365) >= 61 || mod(t,365) <= 151 % Autumn rains
        if rand(1) > 0.95
            r_am = 0.02506 + (-1 + (1+1)*rand(1))*0.01560596; % The downpour event occurs
        end
    elseif mod(t,365) >= 152 || mod(t,365) <= 262 % Winter rains
        if rand(1) > 0.97
            r_am = 0.02016 + (-1 + (1+1)*rand(1))*0.01085544; % The downpour event occurs
        end
    else % Spring rains
        if rand(1) > 0.93
            r_am = 0.02162 + (-1 + (1+1)*rand(1))*0.01228265; % The downpour event occurs
        end
    end
end


end % Function