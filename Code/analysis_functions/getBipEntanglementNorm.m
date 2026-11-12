function [N2mapped] = getBipEntanglementNorm(N2)

    %!!! Calculate the mapped bipartite entanglement norm
    %
    %    input: N2: The raw connected correlation function norm
    %    output: N2mapped: The mapped connected correlation function norm
    
    % Check input parameters
    arguments
        N2 (1,1) double
    end

    % The y which solves the following equation in dependence of the input
    % norm is required to calculate the modified norm
    syms y
    eqn = y^2 - 6*y + (5 - 12*N2) == 0;
    X = solve(eqn, y, 'Real', true);
    
    % Filter solutions
    rangeSol = X(isAlways(X >= -1 & X <= 1));
    sol = double(rangeSol(1));

    % Calculate the norm based on the suitable solution
    N2mapped = (1-sol)/2;
end