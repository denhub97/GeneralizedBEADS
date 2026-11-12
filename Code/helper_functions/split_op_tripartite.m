% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [Op12, Op13, Op23, Op123] = split_op_tripartite(Op)

    %!!! Split a 3-qubit operator into its bipartite and the remaining
    %    components for further analysis of the correlation structure. The
    %    operator is expanded in the Pauli basis, Op = 1/8 * sum_P <P> P
    %    over all Pauli strings P, and the contributions of the Pauli
    %    strings which act non-trivially on exactly two specific qubits are
    %    collected in the bipartite components. The remainder (identity,
    %    single-qubit and three-qubit contributions) is returned as the
    %    fourth output
    %
    %    input: Op: The 8x8 operator to be split, e.g., a density operator
    %               or one of its correlation structures
    %
    %    output: Op12/Op13/Op23: The components of the operator acting
    %                            non-trivially only on qubits 1 and 2, 1
    %                            and 3, or 2 and 3, respectively
    %            Op123: The remaining part of the operator,
    %                       Op - Op12 - Op13 - Op23

    % Check that the operator is 8x8 (3 qubits)
    arguments
        Op (8,8) {mustBeNumeric}
    end

    % Define Paulis
    paulis = {[0 1; 1 0], [0 -1i; 1i 0], [1 0; 0 -1]};
    unit_matrix = eye(2);

    % Component on qubits 1 and 2: sum over all two-qubit Pauli strings
    % with identity on qubit 3
    Op12 = zeros(size(Op));
    for j = 1:3
        for k = 1:3
            Op12 = Op12+ ...
                .125*trace(Op*kron(paulis{j},kron(paulis{k},unit_matrix)))*...
                kron(paulis{j},kron(paulis{k},unit_matrix));
        end
    end

    % Component on qubits 1 and 3: sum over all two-qubit Pauli strings
    % with identity on qubit 2
    Op13 = zeros(size(Op));
    for j = 1:3
        for k = 1:3
            Op13 = Op13+ ...
                .125*trace(Op*kron(paulis{j},kron(unit_matrix,paulis{k})))*...
                kron(paulis{j},kron(unit_matrix,paulis{k}));
        end
    end

    % Component on qubits 2 and 3: sum over all two-qubit Pauli strings
    % with identity on qubit 1
    Op23 = zeros(size(Op));
    for j = 1:3
        for k = 1:3
            Op23 = Op23+ ...
                .125*trace(Op*kron(unit_matrix,kron(paulis{j},paulis{k})))*...
                kron(unit_matrix,kron(paulis{j},paulis{k}));
        end
    end

    % The remainder after subtracting the bipartite components
    Op123 = Op-Op12-Op13-Op23;
end
