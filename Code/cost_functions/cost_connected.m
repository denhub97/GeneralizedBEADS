% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [EntNorm] = cost_connected(DensityOp, Dim)

    %!!! Cost function for calculating the connected correlation function
    %    norm of a 2-qubit or 3-qubit density operator. The density operator
    %    is expanded in the Pauli basis and all contributions which can be
    %    expressed by lower-order correlations (local expectation values and
    %    lower-order connected correlations) are subtracted. The Frobenius
    %    norm of the remainder is the norm of the highest-order connected
    %    correlations. In contrast to cost_connected_tripartite, the Pauli
    %    operators are constructed at runtime
    %
    %    input: DensityOp: The density operator describing the state
    %           Dim: The dimension of the density operator (4 for two
    %                qubits, 8 for three qubits)
    %
    %    output: EntNorm: The connected correlation function norm

    % Pauli matrices and identity
    paulis = {[0 1; 1 0], [0 -1i; 1i 0], [1 0; 0 -1]};
    unit_matrix = eye(2);

    switch Dim

        % For 2 qubit-systems
        case 4

            % Initialize the contributions which are subtracted
            double_contributions = zeros(Dim,Dim);
            single_contributions1 = zeros(Dim,Dim);
            single_contributions2 = zeros(Dim,Dim);

            % Product of the local expectation values of both qubits
            for j = 1:3
                for k = 1:3
                    double_contributions = double_contributions + ...
                        trace(DensityOp*kron(paulis{j},unit_matrix))*...
                        trace(DensityOp*kron(unit_matrix,paulis{k}))*...
                        kron(paulis{j},paulis{k});
                end
            end

            % Local contributions with expectation values of the first qubit
            for m = 1:3
                single_contributions1 = single_contributions1 + ...
                    trace(DensityOp*kron(paulis{m},unit_matrix))*...
                    kron(paulis{m},unit_matrix);
            end

            % Local contributions with expectation values of the second qubit
            for n = 1:3
                single_contributions2 = single_contributions2 + ...
                    trace(DensityOp*kron(unit_matrix,paulis{n}))* ...
                    kron(unit_matrix,paulis{n});
            end

            % Part of the density operator which can be written in terms of
            % local expectation values, including the identity contribution
            rho_nonconnected = 0.25*double_contributions + ...
                               0.25*single_contributions1 + ...
                               0.25*single_contributions2 + ...
                               0.25*trace(DensityOp)*eye(4);

            % The remainder consists of the connected correlations only,
            % compute its Frobenius norm, only consider real part
            rho_tilde = DensityOp - rho_nonconnected;
            EntNorm = real(norm(rho_tilde,'fro'));

        % For 3-qubit systems
        case 8

            % Initialize the contributions which are subtracted
            local_contributions1 = zeros(Dim,Dim);
            local_contributions2 = zeros(Dim,Dim);
            local_contributions3 = zeros(Dim,Dim);
            double_contributions1 = zeros(Dim,Dim);
            double_contributions2 = zeros(Dim,Dim);
            double_contributions3 = zeros(Dim,Dim);
            triple_contributions = zeros(Dim,Dim);
            E12_contribution = zeros(Dim,Dim);
            E13_contribution = zeros(Dim,Dim);
            E23_contribution = zeros(Dim,Dim);

            % Local contributions with expectation values of one single qubit

            % Qubit 1
            for m = 1:3
                local_contributions1 = local_contributions1 + ...
                    trace(DensityOp*kron(kron(paulis{m},unit_matrix),unit_matrix))*...
                    kron(kron(paulis{m},unit_matrix),unit_matrix);
            end

            % Qubit 2
            for n = 1:3
                local_contributions2 = local_contributions2 + ...
                    trace(DensityOp*kron(kron(unit_matrix,paulis{n}),unit_matrix))* ...
                    kron(kron(unit_matrix,paulis{n}),unit_matrix);
            end

            % Qubit 3
            for p = 1:3
                local_contributions3 = local_contributions3 + ...
                    trace(DensityOp*kron(kron(unit_matrix,unit_matrix),paulis{p}))* ...
                    kron(kron(unit_matrix,unit_matrix),paulis{p});
            end

            % Bipartite
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
            % Single-body term products
            for k = 1:3
                for l = 1:3
                    for m = 1:3
                        triple_contributions = triple_contributions + ...
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
                        E12_contribution = E12_contribution + ...
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
                        E13_contribution = E13_contribution + ...
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
                        E23_contribution = E23_contribution + ...
                            (trace(DensityOp*kron(kron(unit_matrix,paulis{l}),paulis{m}))-...
                            (trace(DensityOp*kron(kron(unit_matrix,paulis{l}),unit_matrix))*...
                            trace(DensityOp*kron(kron(unit_matrix,unit_matrix),paulis{m}))))* ...
                            trace(DensityOp*kron(kron(paulis{k},unit_matrix),unit_matrix))*...
                            kron(kron(paulis{k},paulis{l}),paulis{m});
                    end
                end
            end

            % Part of the density operator which can be written in terms of
            % lower-order correlations, including the identity contribution
            rho_nonconnected = 0.125*triple_contributions + ...
                               0.125*double_contributions1 + ...
                               0.125*double_contributions2 + ...
                               0.125*double_contributions3 + ...
                               0.125*local_contributions1 + ...
                               0.125*local_contributions2 + ...
                               0.125*local_contributions3 + ...
                               0.125*E12_contribution + ...
                               0.125*E13_contribution + ...
                               0.125*E23_contribution + ...
                               0.125*trace(DensityOp)*eye(8);

            % The remainder consists of the tripartite connected
            % correlations only, compute its Frobenius norm, only consider
            % real part
            rho_tilde = DensityOp - rho_nonconnected;
            EntNorm = real(norm(rho_tilde,'fro'));

        % Higher qubit systems are not supported
        otherwise
            disp('Connected correlation norm not implemented for higher qubit systems')

    end
end
