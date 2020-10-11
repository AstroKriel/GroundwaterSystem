function J = JacobianInexact(F, u, F0, r)
    % INEXACT JACOBIAN that works 
    N = size(u,1);
    J = zeros(N);
    
    if norm(u,2) == 0
        h = sqrt(eps);
    else
        h = sqrt(eps)*norm(u,2);
    end
    
    fprintf('Progress:\n-----------\n') % loading bar. 
    
    e = zeros(N,1);
    for j = 1:N
        if 1 == mod(j, round(N/10))
            fprintf('*')
        end
        e(j) = 1.0;
        J(:,j) = (F(u+h*e) - F0)/h;
        e(j) = 0.0;
    end

    fprintf('\nDone!\n')
    
    J = sparse(J);
end