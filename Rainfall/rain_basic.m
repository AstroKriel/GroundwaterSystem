function r_am = rain_basic(t, cfg)
% This calculates the rainfall on a given day using the basic cosine model
% r_am refers to rainfall amount
% t refers to the day
% The coefficient of 0.1515 is half the range of 2018's rainfall
% kP, r and N are used for consistency with the other functions used
% Author: Sam

yearRain = cfg.yearRainAverage;

% Calculate the rainfall for a day
if isfield(cfg, 'droughtStart')
    if t >= cfg.droughtStart && t <= cfg.droughtStart + cfg.droughtLength % Enforce drought conditions
        if rand < 0.01 % Allow a small chance of reduced rain
            r_am = 0.5*yearRain*(1 + cos(2*pi*t/365))/365;
        else
            r_am = 0;
        end
    else 
        r_am = yearRain*(1 + cos(2*pi*t/365))/365;
    end
else
    r_am = yearRain*(1 + cos(2*pi*t/365))/365;
end

end