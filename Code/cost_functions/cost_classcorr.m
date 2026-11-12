% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [Cost, P0, P1] = cost_classcorr(UCoeff, DensityOp, TargetQubit, DimSingleQubit, N)

    %!!! Cost function for maximizing the classical correlations of a
    %    bipartite or tripartite state. A single-qubit projective
    %    measurement is performed on the target qubit and the conditional
    %    entropy of the remaining qubits is computed as the cost
    %
    %    input: UCoeff: Coefficients of the Hermitian generator of the
    %                   single-qubit unitary defining the projector
    %                   orientation
    %           DensityOp: The density operator describing the state
    %           TargetQubit: The qubit which is measured
    %           DimSingleQubit: Dimension of the single-qubit Hilbert space
    %           N: Number of qubits of the system (2 or 3)
    %
    %    output: Cost: The conditional entropy after the measurement
    %            P0/P1: The single-qubit projectors on the measured qubit

    % Two-qubit system
    if N == 2

        % Generate 1 qubit unitary matrix
        hermitian = generate_hermitian(UCoeff, DimSingleQubit);
        unitary = expm(1i*hermitian);

        % Rotate the computational basis states to obtain the projectors
        P0 = unitary*[1;0];
        P0 = P0*P0';

        P1 = unitary*[0;1];
        P1 = P1*P1';

        % Extend the projectors to the full system
        if TargetQubit == 1
            P0_fullsys = kron(P0, eye(2));
            P1_fullsys = kron(P1, eye(2));
        elseif TargetQubit == 2
            P0_fullsys = kron(eye(2), P0);
            P1_fullsys = kron(eye(2), P1);
        end

        % Compute the probabilities and the post-measurement states of the
        % unmeasured qubit for both outcomes
        p0 = trace(P0_fullsys*DensityOp);
        outcome0 = partial_trace(P0_fullsys*DensityOp*P0_fullsys, TargetQubit)/p0;
        p1 = trace(P1_fullsys*DensityOp);
        outcome1 = partial_trace(P1_fullsys*DensityOp*P1_fullsys, TargetQubit)/p1;

        % Calculate the conditional entropy in terms of the outcome
        % eigenvalue decomposition, eigenvalues smaller than tolerance do
        % not contribute to the entropy (evaluate to 0). This circumvents
        % the problem of logm() = NaN, when rho does not have full rank.

        % The numerical tolerance for eigenvalues
        tolerance = 1e-8;

        % Compute the eigendecomposition of outcome corresp. to P0
        [~,D] = eig(outcome0);
        eigvalues_0 = diag(D);

        % Compute the von Neumann entropy of outcome corresp. to P0
        s0 = 0;
        for i = 1:size(eigvalues_0,1)
            % Continue if eigenvalue is smaller than tolerance
            if abs(eigvalues_0(i)) < tolerance
                continue;
            end

            s0 = s0 + eigvalues_0(i)*log(eigvalues_0(i));
        end

        % Weight with the outcome probability
        s0 = p0*(-s0);

        % Compute the eigendecomposition of outcome corresp. to P1
        [~,D] = eig(outcome1);
        eigvalues_1 = diag(D);

        % Compute the von Neumann entropy of outcome corresp. to P1
        s1 = 0;
        for i = 1:size(eigvalues_1,1)
            % Continue if eigenvalue is smaller than tolerance
            if abs(eigvalues_1(i)) < tolerance
                continue;
            end

            s1 = s1 + eigvalues_1(i)*log(eigvalues_1(i));
        end

        % Weight with the outcome probability
        s1 = p1*(-s1);

        % Conditional entropy, only consider real part
        Cost = real(s0 + s1);

    % Three-qubit system
    elseif N == 3

        % Generate 1 qubit unitary matrix
        hermitian = generate_hermitian(UCoeff, DimSingleQubit);
        unitary = expm(1i*hermitian);

        % Rotate the computational basis states to obtain the projectors
        P0 = unitary*[1;0];
        P0 = P0*P0';
        P1 = unitary*[0;1];
        P1 = P1*P1';

        % Extend the projectors to the full system
        if TargetQubit == 1
            P0_fullsys = kron(P0, eye(4));
            P1_fullsys = kron(P1, eye(4));
        elseif TargetQubit == 2
            P0_fullsys = kron(eye(2), kron(P0, eye(2)));
            P1_fullsys = kron(eye(2), kron(P1, eye(2)));
        elseif TargetQubit == 3
            P0_fullsys = kron(eye(4), P0);
            P1_fullsys = kron(eye(4), P1);
        end

        % Compute the probabilities and the post-measurement states of the
        % unmeasured qubits for both outcomes
        p0 = trace(P0_fullsys*DensityOp);
        outcome0 = partial_trace(P0_fullsys*DensityOp*P0_fullsys, TargetQubit)/p0;
        p1 = trace(P1_fullsys*DensityOp);
        outcome1 = partial_trace(P1_fullsys*DensityOp*P1_fullsys, TargetQubit)/p1;

        % Calcuate the conditional entropy in terms of the outcome
        % eigenvalue decomposition, eigenvalues smaller than tolerance do
        % not contribute to the entropy (evaluate to 0). This circumvents
        % the problem of logm() = NaN, when rho does not have full rank.

        % The numerical tolerance for eigenvalues
        tolerance = 1e-8;

        % Compute the eigendecomposition of outcome corresp. to P0
        [~,D] = eig(outcome0);
        eigvalues_0 = diag(D);

        % Compute the von Neumann entropy of outcome corresp. to P0
        s0 = 0;
        for i = 1:size(eigvalues_0,1)
            % Continue if eigenvalue is smaller than tolerance
            if abs(eigvalues_0(i)) < tolerance
                continue;
            end

            s0 = s0 + eigvalues_0(i)*log(eigvalues_0(i));
        end

        % Weight with the outcome probability
        s0 = p0*(-s0);

        % Compute the eigendecomposition of outcome corresp. to P1
        [~,D] = eig(outcome1);
        eigvalues_1 = diag(D);

        % Compute the von Neumann entropy of outcome corresp. to P1
        s1 = 0;
        for i = 1:size(eigvalues_1,1)
            % Continue if eigenvalue is smaller than tolerance
            if abs(eigvalues_1(i)) < tolerance
                continue;
            end

            s1 = s1 + eigvalues_1(i)*log(eigvalues_1(i));
        end

        % Weight with the outcome probability
        s1 = p1*(-s1);

        % Conditional entropy, only consider real part
        Cost = real(s0 + s1);
    end
end
