function r_av = average_rain(cfg, t)
% This function calculates the average yearly rainfall for the coming year
% using the least squares cosine model at the start of each year
% r_av refers to the average rainfall amount
% d refers to the day
% kP, r and N are used for consistency with the other functions used

% Determine the year
y = floor(t/365);

% Set the days that rainfall will be predicted for
days = 1 + (y-1)*365:1:365 + (y-1)*365;

% Calculate the predicted values
r_pred = 1.883*10^(-3) + 2.702*10^(-5)*cos(2*pi*days/1825) + 1.023*10^(-3)...
    *cos(2*pi*days/365) - 1.719*10^(-4)*cos(2*pi*days/90) - 6.189*10^(-5)*cos(2*pi*days/30);

% Calculate the average rainfall
r_av = mean(r_pred);
