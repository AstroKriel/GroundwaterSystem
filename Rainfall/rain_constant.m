function r_am = rain_constant(~, cfg)
% RAIN_CONSTANT takes inputs t which is the time given in days, where r_d
% is the rainfall amount. it is assumed that this function is only called
% on each new day. 

% Average yearly rainfall 
R = cfg.yearRainAverage;

r_am = R/365;

end