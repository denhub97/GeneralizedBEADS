function [Concurrence] = get_concurrence(DensityOp)

    %!!! Calculate the concurrence of a bipartite density operator
    %   
    %    input: DensityOp: The density operator describing the state
    %
    %    output: Concurrence: The concurrence

    % Check input arguments
    arguments
        DensityOp (4,4) {mustBeNumeric}
    end

    % Pauli-y definition
    pauli_y = [0 -1i; 1i 0];
    
    % Spin-flipped density matrix
    rho_tilde = kron(pauli_y,pauli_y) * conj(DensityOp) * kron(pauli_y,pauli_y);
    
    % Non-hermitian matrix of which we compute sqrt(eigenvalues)
    non_hermitian = DensityOp * rho_tilde;
    
    % Compute the eigendecomposition
    [~, eigenvalues] = eig(non_hermitian);
    lambda = sqrt(diag(eigenvalues));
    
    % Sort the eigenvalues in descending order
    lambda = sort(lambda, 'descend');
    
    % Choose max value 0 or the difference of sqrt(eigenvalues) for
    % concurrence
    Concurrence = max(0, real(lambda(1) - lambda(2) - lambda(3) - lambda(4)));
end