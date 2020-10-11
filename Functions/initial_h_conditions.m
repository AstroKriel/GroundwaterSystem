function h = initial_h_conditions(h_bot, h_top, Z, L2)
    % Using the initial conditions given by hbot and htop we have a gradient of
    % flow.
    
    % Mitch: this function isn't very complex and only seems to be called once,
    % so it might be worth absorbing into the code that is calling it.

    h = h_bot .* ones(size(Z)) + ((h_top - h_bot)/L2) .* Z;
    h = h(:);
end