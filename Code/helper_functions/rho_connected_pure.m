% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [RhoConn, RhoComp] = rho_connected_pure(DensityOp, Dim)

    %!!! Split the density operator of a pure state of a 2-qubit or 3-qubit
    %    system into its connected and its compound part. The density
    %    operator is expanded in the Pauli basis and all contributions
    %    which can be expressed as products of local expectation values
    %    are subtracted: the unconnected product contributions of two or
    %    three qubits and, for three qubits, the contributions composed of
    %    a two-qubit connected correlation and the expectation value of the
    %    third qubit. The connected part (for visualization) retains the
    %    identity, the local (single-qubit) contributions, and the
    %    connected correlations. The compound part is the difference to the
    %    density operator, i.e., the subtracted contributions. For a mixed
    %    state, the results of the individual pure states of a
    %    decomposition are summed with their weights
    %
    %    input: DensityOp: The density operator describing the state (4x4
    %                      or 8x8)
    %           Dim: The dimension of the density operator (4 for two
    %                qubits, 8 for three qubits)
    %
    %    output: RhoConn: The connected part of the density operator
    %            RhoComp: The compound part of the density operator,
    %                     DensityOp - RhoConn

    % Check input parameters, the density operator must be a numerical
    % 4x4 or 8x8 matrix with matching dimension
    arguments
        DensityOp (:,:) {mustBeNumeric}
        Dim (1,1) {mustBeNumeric}
    end
    assert(size(DensityOp,1) == size(DensityOp,2) && ismember(size(DensityOp,1), [4 8]), ...
           "rho_connected_pure:invalidDimension", ...
           "DensityOp must be a square matrix of size 4 or 8, but is %d x %d.", ...
           size(DensityOp,1), size(DensityOp,2));
    assert(Dim == size(DensityOp,1), ...
           "rho_connected_pure:dimensionMismatch", ...
           "Dim (%d) does not match the size of DensityOp (%d).", ...
           Dim, size(DensityOp,1));

    % Pauli matrices and identity
    paulis = {[0 1; 1 0], [0 -1i; 1i 0], [1 0; 0 -1]};
    unit_matrix = eye(2);

    switch Dim

        % For 2 qubit-systems
        case 4

            % Initialize the unconnected contribution of both qubits
            double_contributions = zeros(Dim,Dim);

            % Product of the local expectation values of both qubits
            for j = 1:3
                for k = 1:3
                    double_contributions = double_contributions + ...
                        trace(DensityOp*kron(paulis{j},unit_matrix))*...
                        trace(DensityOp*kron(unit_matrix,paulis{k}))*...
                        kron(paulis{j},paulis{k});
                end
            end

            % Subtract the unconnected contribution, the remainder is the
            % connected part
            RhoConn = DensityOp - 0.25*double_contributions;
            RhoComp = DensityOp-RhoConn;

        % For 3-qubit systems
        case 8

            % Initialize the contributions which are subtracted
            double_contributions1 = zeros(Dim,Dim);
            double_contributions2 = zeros(Dim,Dim);
            double_contributions3 = zeros(Dim,Dim);
            triple_contributions = zeros(Dim,Dim);
            E12_contribution = zeros(Dim,Dim);
            E13_contribution = zeros(Dim,Dim);
            E23_contribution = zeros(Dim,Dim);

            % Bipartite, products of the local expectation values of two qubits
            % Qubits 1 and 2
            for j = 1:3
                for k = 1:3
                    double_contributions1 = double_contributions1 + ...
                        trace(DensityOp*kron(kron(paulis{j},unit_matrix),unit_matrix))*...
                        trace(DensityOp*kron(kron(unit_matrix,paulis{k}),unit_matrix))*...
                        kron(kron(paulis{j},paulis{k}),unit_matrix);
                end
            end

            % Qubits 1 and 3
            for f = 1:3
                for g = 1:3
                    double_contributions2 = double_contributions2 + ...
                        trace(DensityOp*kron(kron(paulis{f},unit_matrix),unit_matrix))*...
                        trace(DensityOp*kron(kron(unit_matrix,unit_matrix),paulis{g}))*...
                        kron(kron(paulis{f},unit_matrix),paulis{g});
                end
            end

            % Qubits 2 and 3
            for h = 1:3
                for i = 1:3
                    double_contributions3 = double_contributions3 + ...
                        trace(DensityOp*kron(kron(unit_matrix,paulis{h}),unit_matrix))*...
                        trace(DensityOp*kron(kron(unit_matrix,unit_matrix),paulis{i}))*...
                        kron(kron(unit_matrix,paulis{h}),paulis{i});
                end
            end
    
            % Tripartite
            % Product of the local expectation values of all three qubits
            for k = 1:3
                for l = 1:3
                    for m = 1:3
                        triple_contributions = triple_contributions+...
                            trace(DensityOp*kron(kron(paulis{k},unit_matrix),unit_matrix))*...
                            trace(DensityOp*kron(kron(unit_matrix,paulis{l}),unit_matrix))*...
                            trace(DensityOp*kron(kron(unit_matrix,unit_matrix),paulis{m}))* ...
                            kron(kron(paulis{k},paulis{l}),paulis{m});
                    end
                end
            end

            % Contributions from connected correlations E_kl * expectation
            % value of a single qubit m, where the connected two-qubit
            % correlation is E_kl = <P_k P_l> - <P_k><P_l>
            % E12 x Q3
            for k = 1:3
                for l = 1:3
                    for m = 1:3
                        E12_contribution = E12_contribution+...
                            (trace(DensityOp*kron(kron(paulis{k},paulis{l}),unit_matrix))-...
                            (trace(DensityOp*kron(kron(paulis{k},unit_matrix),unit_matrix))*...
                            trace(DensityOp*kron(kron(unit_matrix,paulis{l}),unit_matrix))))* ...
                            trace(DensityOp*kron(kron(unit_matrix,unit_matrix),paulis{m}))*...
                            kron(kron(paulis{k},paulis{l}),paulis{m});
                    end
                end
            end

            % E13 x Q2
            for k = 1:3
                for l = 1:3
                    for m = 1:3
                        E13_contribution = E13_contribution+...
                            (trace(DensityOp*kron(kron(paulis{k},unit_matrix),paulis{m}))-...
                            (trace(DensityOp*kron(kron(paulis{k},unit_matrix),unit_matrix))*...
                            trace(DensityOp*kron(kron(unit_matrix,unit_matrix),paulis{m}))))* ...
                            trace(DensityOp*kron(kron(unit_matrix,paulis{l}),unit_matrix))*...
                            kron(kron(paulis{k},paulis{l}),paulis{m});
                    end
                end
            end

            % E23 x Q1
            for k = 1:3
                for l = 1:3
                    for m = 1:3
                        E23_contribution = E23_contribution+...
                            (trace(DensityOp*kron(kron(unit_matrix,paulis{l}),paulis{m}))-...
                            (trace(DensityOp*kron(kron(unit_matrix,paulis{l}),unit_matrix))*...
                            trace(DensityOp*kron(kron(unit_matrix,unit_matrix),paulis{m}))))* ...
                            trace(DensityOp*kron(kron(paulis{k},unit_matrix),unit_matrix))*...
                            kron(kron(paulis{k},paulis{l}),paulis{m});
                    end
                end
            end

            % Subtract all compound contributions, the
            % remainder is the connected part
            RhoConn = DensityOp-(0.125*triple_contributions+...
                                 0.125*double_contributions1 + ...
                                 0.125*double_contributions2 + ...
                                 0.125*double_contributions3+...
                                 0.125*E12_contribution+...
                                 0.125*E13_contribution+...
                                 0.125*E23_contribution);

            % The compound part is the difference to the density operator
            RhoComp = DensityOp-RhoConn;
    end
end
