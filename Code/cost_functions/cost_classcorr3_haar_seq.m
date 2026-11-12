% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [Cost, P0FullSys, P1FullSys, P0] = cost_classcorr3_haar_seq(UCoeff, DensityOp, TargetQubit)

    %!!! Cost function for maximizing the classical correlations of a
    %    tripartite state. A single-qubit projective measurement is
    %    performed on the target qubit and the conditional entropy of the
    %    remaining qubits is computed as the cost
    %
    %    input: UCoeff: Euler angles [theta, phi, lambda] defining the
    %                   single-qubit unitary which sets the projector
    %                   orientation
    %           DensityOp: The density operator describing the tripartite
    %                      state
    %           TargetQubit: The qubit which is measured
    %
    %    output: Cost: The conditional entropy after the measurement
    %            P0FullSys/P1FullSys: The projectors extended to the full
    %                                 3-qubit system
    %            P0: The single-qubit projector corresponding to outcome 0

    % Generate single-qubit unitary matrix
    unitary1 = [cos(UCoeff(1)/2), -exp(1i*UCoeff(3)) * sin(UCoeff(1)/2);
                exp(1i*UCoeff(2)) * sin(UCoeff(1)/2), exp(1i*(UCoeff(2) + UCoeff(3))) * cos(UCoeff(1)/2)];

    % Rotate the computational basis states to obtain the projectors
    P0 = unitary1*[1;0];
    P0 = P0*P0';
    P1 = unitary1*[0;1];
    P1 = P1*P1';

    % Extend the projectors to the full system
    if TargetQubit == 1
        P0FullSys = kron(P0, eye(4));
        P1FullSys = kron(P1, eye(4));
    elseif TargetQubit == 2
        P0FullSys = kron(eye(2), kron(P0, eye(2)));
        P1FullSys = kron(eye(2), kron(P1, eye(2)));
    elseif TargetQubit == 3
        P0FullSys = kron(eye(4), P0);
        P1FullSys = kron(eye(4), P1);
    end

    % Compute the probabilities and the post-measurement states of the
    % unmeasured qubits for both outcomes
    p0 = trace(P0FullSys*DensityOp);
    outcome0 = partial_trace(P0FullSys*DensityOp*P0FullSys, TargetQubit) / p0;

    p1 = trace(P1FullSys*DensityOp);
    outcome1 = partial_trace(P1FullSys*DensityOp*P1FullSys, TargetQubit) / p1;

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

        s0 = s0 + eigvalues_0(i) * log(eigvalues_0(i));
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

        s1 = s1 + eigvalues_1(i) * log(eigvalues_1(i));
    end

    % Weight with the outcome probability
    s1 = p1*(-s1);

    % Conditional entropy, only consider real part
    Cost = real(s0 + s1);
end
