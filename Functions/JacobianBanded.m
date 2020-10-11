function J = JacobianBanded(Ffunc, x, Fx0, r)
    % Banded Jacobian.
 	% Author: Neco Kriel
 	% Input:
 	%   Ffunc     Function that calculates flow
 	%   x         previous iteration's position
 	%   Fx0       previous solution
 	%   r         Number of nodes in a column for our domain
 	% Output: Jacobian solution
    
    bw = 2*r + 1;          % bandwidth of our system
    N  = length(x);        % number of nodes. old version: size(x,1)
    J  = spalloc(N,N,5*N); % initialise Jacobian

    % spalloc pre allocates memory assuming an NxN matrix with atmost 5 x N
    % elements which is true for us (each column has max of 5 values. This
    % reduced time by 10%.

    % Choose the size of the shift parameter.
    ynorm = norm(x,2);
    if ynorm == 0
        h = sqrt(eps);
    else
        h = ynorm*sqrt(eps);
    end
    
    % Calculating Jacobian column by column.
    s = zeros(N,1);
    for j = 1:bw
        
        s(j:bw:N) = 1;                % grab column
        h_norm = h / norm(s,2);       % norm step
        Fpert  = Ffunc(x + h_norm*s); % perturb soln
        J_col  = (Fpert-Fx0)/h_norm;  % calc jacobian column
        s(j:bw:N) = 0;                % reset column grabbing vector
        
        col       = ( j:bw:N )';
        start_row = ( max(col-r,1) )';
        end_row   = ( min(col+r,N) )';
        
        for i = 1:length(start_row)
            J(start_row(i):end_row(i), col(i)) = J_col(start_row(i):end_row(i));
        end
    end
end
