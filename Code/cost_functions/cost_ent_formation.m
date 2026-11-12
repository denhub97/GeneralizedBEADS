function [EoF, States] = cost_ent_formation(UCoeff, Eigenstates, RankRho, Dim)

    %!!! Cost function for minimizing the average entanglement of a pure
    %    state decomposition of a bipartite state. A unitary mixes the
    %    eigenstates of the density operator to generate a new pure state
    %    ensemble and the average entanglement entropy of the ensemble is
    %    computed as the cost
    %
    %    input: UCoeff: Coefficients of the Hermitian generator of the
    %                   unitary mixing the eigenstates
    %           Eigenstates: The weighted eigenstates of the density
    %                        operator, stored as columns of the matrix
    %           RankRho: The number of states in the decomposition
    %           Dim: Dimension of the Hilbert space
    %
    %    output: EoF: The average entanglement of the pure state ensemble
    %            States: The cell array of pure states (including their
    %                    weights) constituting the new decomposition

    % Generate the unitary matrix mixing the eigenstates
    hermitian = generate_hermitian(UCoeff, RankRho);
    unitary = expm(1i*hermitian);

    % With the eigendecomposition as input, the eigenstates are the columns
    % of the Eigenstates matrix. Construct the new states as linear
    % combinations of the eigenstates
    States = cell(RankRho, 1);
    for k = 1:RankRho
        phi_k = zeros(Dim, 1);
        for l = 1:RankRho
            phi_k = phi_k + unitary(k,l)*Eigenstates(:,l);
        end
        States{k,1} = phi_k;
    end

    % Initialize the state weights and entanglement values
    list_of_norms = zeros(1, RankRho);
    list_ent_of_formation = zeros(1, RankRho);

    % Compute the entanglement of each pure state
    for index = 1:RankRho
        % The squared norm of the state is its weight in the ensemble
        list_of_norms(index) = norm(States{index});

        % Normalize the state and compute the reduced density operator of
        % the first qubit
        norm_vec = States{index}/norm(States{index});
        rho_temp = norm_vec*norm_vec';
        rho_reduced = partial_trace(rho_temp, 1);

        % Turn off the warning about negative eigenvalues for logm because
        % it is handeled later on in the code. Remember the previous
        % warning state to restore it afterwards
        warningState = warning('off','all');

        % Entanglement entropy of the pure state
        ent_formation = (-1)*trace(rho_reduced*logm(rho_reduced));

        % Ensure that the entropy is numerical
        if isnan(ent_formation)
            ent_formation = 0;
        end

        list_ent_of_formation(index) = ent_formation;
    end

    % Restore the previous warning state
    warning(warningState);

    % Average the entanglement of the states weighted by their probability,
    % only consider real part
    EoF = real(sum((list_of_norms.^2) .* list_ent_of_formation));
end
