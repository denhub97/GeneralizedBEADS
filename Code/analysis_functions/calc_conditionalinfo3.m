function [CondInfo] = calc_conditionalinfo3(DensityOp, P00, P01, P10, P11, TargetQubits)

    %!!! Function for calculating the conditional entropy for a joint local
    %    measurement of two qubits in a 3-qubit system
    %
    %    input:  DensityOp: The density operator describing the tripartite state
    %            P00/P01/P10/P11: Full set of projectors, must be 8x8
    %                            matrices
    %            TargetQubits: The list of target qubits being measured
    % 
    %    output: CondInfo: Conditional entropy

    % Check inputs, ensure correct data types
    arguments
        DensityOp (8,8) {mustBeNumeric}
        P00 (8,8) {mustBeNumeric}
        P01 (8,8) {mustBeNumeric}
        P10 (8,8) {mustBeNumeric}
        P11 (8,8) {mustBeNumeric}
        TargetQubits (1,2) {mustBeNumeric}
    end

    % Sort the list of target qubits in descending order for efficient
    % partial trace operations
    TargetQubits = sort(TargetQubits,'descend');

    % Collect projectors, probabilities and outcomes
    projectors = {P00, P01, P10, P11};
    probabilities = zeros(1,4);
    outcomes = cell(1,4);

    % Compute the probability and outcome for each projector
    for pp = 1:numel(probabilities)
        % The current projector
        P = projectors{pp};
        
        probabilities(pp) = trace(P * DensityOp);
        outcomes{pp} = partial_trace(partial_trace(P * DensityOp * P, TargetQubits(1)), TargetQubits(2)) / probabilities(pp);
    end
    
    % Calcuate the conditional entropy in terms of the outcome
    % eigenvalue decomposition, eigenvalues smaller than tolerance do not
    % contribute to the entropy (evaluate to 0). This circumvents the
    % problem of logm() = NaN, when rho does not have full rank.
    
    % The numerical tolerance for eigenvalues
    tolerance = 1e-8;
    
    CondInfo = 0;
    for ii = 1:numel(probabilities)
        p = probabilities(ii);
        outcome = outcomes{ii};

        % Compute the eigendecomposition of outcome corresp. to P11
        [~,D] = eig(outcome);
        eigvalues = diag(D);
    
        % Compute the von Neumann entropy of outcome corresp. to P11
        S = 0;
        for jj = 1:size(eigvalues, 1)
            % Continue if eigenvalue is smaller than tolerance
            if abs(eigvalues(jj)) < tolerance
               continue;
            end
        
            S = S + eigvalues(jj) * log(eigvalues(jj));
        
        end
        S = p * (-S);
        
        % Add to conditional entropy, only consider real part
        CondInfo = CondInfo +  real(S);
    end
end