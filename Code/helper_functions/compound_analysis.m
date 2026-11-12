function [RhoSigma] = compound_analysis(DensityOp, RhoComp, RhoConn, Dim)

    %!!! Compute the genuinely correlation-based separable structure of a
    %    density operator. Starting from the compound structure (the
    %    weighted sum of the compound structures of the pure states which
    %    contribute to a mixed quantum state), all contributions which are
    %    uncorrelated in the information-theoretical sense are subtracted,
    %    i.e., contributions which factorize into products of
    %    expectation values of local measurements (local contributions and
    %    unconnected product contributions of two or three qubits). For
    %    three qubits, also the contributions composed of a two-qubit
    %    correlation and the expectation value of the third qubit are
    %    removed. What remains is the separable structure which is
    %    solely due to correlations between the subsystems
    %
    %    input: DensityOp: The density operator describing the state
    %           RhoComp: The compound structure of the density operator
    %                    (weighted sum of the compound structures of the
    %                    contributing pure states)
    %           RhoConn: The connected correlation structure of the
    %                    density operator (weighted sum of the connected
    %                    structures of the contributing pure states)
    %           Dim: The dimension of the density operator (4 for two
    %                qubits, 8 for three qubits)
    %
    %    output: RhoSigma: The genuinely correlation-based separable
    %                      structure of the density operator

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

            % Subtract the unconnected contribution from the compound
            % structure
            RhoSigma = RhoComp - 0.25*double_contributions;

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

            % First subtract the bipartite redundant information
            RhoSigma_bipartite = RhoComp-0.125*double_contributions1 + ...
                                 -0.125*double_contributions2 + ...
                                 -0.125*double_contributions3;

            % Tripartite
            % Single-qubit term products, Qubits 1, 2, and 3
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

            % Contributions from two-qubit correlations E_kl (connected and
            % separable bipartite structure) times the expectation value of
            % the remaining single qubit m
            % E12 x Q3
            for k = 1:3
                for l = 1:3
                    for m = 1:3
                        E12_contribution = E12_contribution+...
                            (trace(RhoConn*kron(kron(paulis{k},paulis{l}),unit_matrix))*...
                            trace(DensityOp*kron(kron(unit_matrix,unit_matrix),paulis{m}))+...
                            ...
                            trace(RhoSigma_bipartite*kron(kron(paulis{k},paulis{l}),unit_matrix))*...
                            trace(DensityOp*kron(kron(unit_matrix,unit_matrix),paulis{m})))*...
                            kron(kron(paulis{k},paulis{l}),paulis{m});
                    end
                end
            end

            % E13 x Q2
            for k = 1:3
                for l = 1:3
                    for m = 1:3
                        E13_contribution = E13_contribution+...
                            (trace(RhoConn*kron(kron(paulis{k},unit_matrix),paulis{m}))*...
                            trace(DensityOp*kron(kron(unit_matrix,paulis{l}),unit_matrix))+...
                            ...
                            trace(RhoSigma_bipartite*kron(kron(paulis{k},unit_matrix),paulis{m}))*...
                            trace(DensityOp*kron(kron(unit_matrix,paulis{l}),unit_matrix)))*...
                            kron(kron(paulis{k},paulis{l}),paulis{m});
                    end
                end
            end

            % E23 x Q1
            for k = 1:3
                for l = 1:3
                    for m = 1:3
                        E23_contribution = E23_contribution+...
                            (trace(RhoConn*kron(kron(unit_matrix,paulis{l}),paulis{m}))*...
                            trace(DensityOp*kron(kron(paulis{k},unit_matrix),unit_matrix))+...
                            ...
                            trace(RhoSigma_bipartite*kron(kron(unit_matrix,paulis{l}),paulis{m}))*...
                            trace(DensityOp*kron(kron(paulis{k},unit_matrix),unit_matrix)))*...
                            kron(kron(paulis{k},paulis{l}),paulis{m});
                    end
                end
            end

            % Subtract the remaining contributions from the structure
            % after the bipartite subtraction
            RhoSigma = RhoSigma_bipartite ...
                       -0.125*triple_contributions ...
                       -0.125*E12_contribution ...
                       -0.125*E13_contribution ...
                       -0.125*E23_contribution;
    end
end
