function [ModDensityOp] = removeIdLinear(DensityOp, RemoveId)

    %!!! Remove the linear operator contributions, i.e., all terms which
    %    only act non-trivially on a single qubit, and optionally also the
    %    identity contribution from the Pauli expansion of the density
    %    operator of a 2-qubit or 3-qubit system. The density operator is
    %    expanded as rho = 1/dim * sum_P <P> P over all Pauli strings P,
    %    and the identity and single-qubit Pauli strings are subtracted
    %
    %    input: DensityOp: The density operator describing the state (4x4
    %                      or 8x8)
    %           RemoveId: Shall the identity contribution be removed as
    %                     well (boolean, default true)
    %
    %    output: ModDensityOp: The density operator without the linear
    %                          contributions (and without the identity
    %                          contribution if requested)

    % Check input parameters
    arguments
        DensityOp (:,:) {mustBeNumeric}
        RemoveId (1,1) logical = true
    end

    % Define Paulis
    paulis = {[0 1; 1 0], [0 -1i; 1i 0], [1 0; 0 -1]};
    unit_matrix = eye(2);
    ModDensityOp = DensityOp;

    % Two-qubit system
    if size(DensityOp,1) == 4

        % Subtract the linear contributions of both qubits
        for j = 1:3
            ModDensityOp = ModDensityOp - 0.25 * ( ...
                           trace(DensityOp * kron(paulis{j}, unit_matrix)) * ...
                                       kron(paulis{j}, unit_matrix)  + ...
                           trace(DensityOp * kron(unit_matrix, paulis{j})) * ...
                                       kron(unit_matrix, paulis{j}) ...
                           );
        end

        % Subtract the identity contribution
        if RemoveId
            ModDensityOp = ModDensityOp - 0.25 * trace(DensityOp * kron(unit_matrix, unit_matrix)) * ...
                                                       kron(unit_matrix, unit_matrix);
        end

    % Three-qubit system
    elseif size(DensityOp,1) == 8

        % Subtract the linear contributions of all three qubits
        for j = 1:3
            ModDensityOp = ModDensityOp - 0.125 * ( ...
                           trace(DensityOp * kron(paulis{j}, kron(unit_matrix, unit_matrix))) * ...
                                       kron(paulis{j}, kron(unit_matrix, unit_matrix))  + ...
                           trace(DensityOp * kron(unit_matrix, kron(paulis{j}, unit_matrix))) * ...
                                       kron(unit_matrix, kron(paulis{j}, unit_matrix))  + ...
                           trace(DensityOp * kron(unit_matrix, kron(unit_matrix, paulis{j}))) * ...
                                       kron(unit_matrix, kron(unit_matrix, paulis{j}))    ...
                           );
        end

        % Subtract the identity contribution
        if RemoveId
            ModDensityOp = ModDensityOp - 0.125 * trace(DensityOp * kron(unit_matrix, kron(unit_matrix, unit_matrix))) * ...
                                                        kron(unit_matrix, kron(unit_matrix, unit_matrix));
        end
    end
end
