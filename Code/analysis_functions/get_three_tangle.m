function [TangleSum, Tangle12, Tangle13, Tangle23, Tangle123] = get_three_tangle(DensityOp)

    %!!! Calculate different tangles of a tripartite pure state
    %   
    %    input: DensityOp: The density operator describing the state, must
    %                      be 8x8 and rank 1 (pure state)
    %
    %    output: TangleSum: The sum of the average 2-tangle and the
    %                       3-tangle
    %            Tangle12/13/23: The 2-tangles of the denoted bipartite
    %                            marginals
    %            Tangle123: The 3-tangle
    
    % Check input arguments
    arguments
        DensityOp (8,8) {mustBeNumeric}
    end

    % Check if the state is really pure
    if rank(DensityOp) ~= 1
        warning("The input state is not pure! The returned three-tangle is potentially incorrect.")
    end

    % Get the reduced density matrices of the bipartite marginals
    rho12 = partial_trace(DensityOp,3);
    rho13 = partial_trace(DensityOp,2);
    rho23 = partial_trace(DensityOp,1);

    % Compute the bipartite concurrences
    C12 = get_concurrence(rho12);
    C13 = get_concurrence(rho13);
    C23 = get_concurrence(rho23);

    % Compute the 2-tangles
    Tangle12 = C12^2;
    Tangle13 = C13^2;
    Tangle23 = C23^2;

    % Compute the average 2-tangle    
    AverageTangle = (Tangle12 + Tangle13 + Tangle23)/3;

    % Get the single-qubit reduced density operator of qubit 1
    rho1 = partial_trace(rho12,2);

    % Compute the concurrence across the bipartition 1:23
    C1_23 = sqrt(2 - 2*real(trace(rho1^2)));

    % Compute the 3-tangle. Follows directly from Coffman–Kundu–Wootters 
    % (CKW) monogamy relation
    Tangle123 = C1_23^2 - (C12^2 + C13^2);

    % Tangle sum
    TangleSum = AverageTangle + Tangle123;
end


