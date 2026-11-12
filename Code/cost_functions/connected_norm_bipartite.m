function [EntNorm, States] = connected_norm_bipartite(UCoeff, Eigenstates, RankRho, Dim)

    %!!! Cost function for minimizing the connected correlation function norm
    %    of a pure state decomposition of a bipartite state. A unitary mixes
    %    the eigenstates of the density operator to generate a new pure
    %    state ensemble and the average connected correlation function norm
    %    of the ensemble is computed as the cost
    %
    %    input: UCoeff: Coefficients of the Hermitian generator of the
    %                   unitary mixing the eigenstates
    %           Eigenstates: The weighted eigenstates of the density
    %                        operator, stored as columns of the matrix
    %           RankRho: The number of states in the decomposition
    %           Dim: Dimension of the Hilbert space
    %
    %    output: EntNorm: The average connected correlation function norm
    %                     of the pure state ensemble
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

    % Initialize the state weights and connected correlation norms
    list_of_norms = zeros(1, RankRho);
    list_ent_rho_connected = zeros(1, RankRho);

    % Compute the connected correlation norm of each pure state
    for index = 1:RankRho
        % The squared norm of the state is its weight in the ensemble
        list_of_norms(index) = norm(States{index});

        % Normalize the state and compute the corresponding density operator
        norm_vec = States{index}/norm(States{index});
        rho_temp = norm_vec*norm_vec';

        % Connected correlation function norm of the pure state
        list_ent_rho_connected(index) = cost_connected(rho_temp, Dim);
    end

    % Average the connected correlation norms of the states weighted by
    % their probability
    EntNorm = sum((list_of_norms.^2).* list_ent_rho_connected);
end
