% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [N3mapped] = getEntanglementNorm(N3)

    %!!! Calculate the mapped tripartite entanglement norm which is sensitive to
    %    distinguishing GHZ- and W- type entanglement
    %    The mapped norm coincides with the three tangle for structurally ideal 
    %    GHZ-entangled states but is smaller or larger if the structure
    %    deviates
    %
    %    input: N3: The original tripartite connected correlation function norm
    %    output: N3mapped: The mapped norm
    
    % Check input parameters
    arguments
        N3 (1,1) double
    end

    % The y which solves the following equation in dependence of the input
    % norm is required to calculate the modified norm
    syms y
    eqn = y^3 - y^2 - 5*y + (5 - 8*N3) == 0;
    X = solve(eqn, y, 'Real', true);

    % Filter solutions
    rangeSol = X(isAlways(X >= -1 & X <= 1));
    sol = double(rangeSol(1));

    % Calculate the norm based on a suitable solution
    N3mapped = (1-sol)/2;
end