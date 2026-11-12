function [Cost, RhoSigmaCl, P0_, P1_] = cost_find_projectors(Angles, NormalVector, P0, P1, RhoSigma, TargetQubit, Qubits, OptStructure)

    %!!! Cost function for finding the single-qubit projectors on the target
    %    qubit which maximize the overlap between the classical correlation
    %    structure accessed by the projective measurement and the separable
    %    structure of a bi- or tripartite density operator. The initial
    %    projectors are rotated about a fixed axis by the angle to be
    %    optimized
    %
    %    input: Angles: Rotation angle of the projectors
    %           NormalVector: Rotation axis (does not need to be normalized)
    %           P0/P1: Initial projectors, extended to the full system
    %           RhoSigma: Separable structure of density operator
    %           TargetQubit: The qubit which is measured
    %           Qubits: Number of qubits of the system (2 or 3)
    %           OptStructure: Structure which the projectors are optimized
    %                         for (integer): 1 = separable structure (the
    %                         overlap is maximized), -1 = entangling
    %                         structure (the overlap is minimized)
    %
    %    output: Cost: Negative overlap shifted by 1 (the overlap is
    %                  multiplied by OptStructure)
    %            RhoSigmaCl: The classical structure of RhoSigma with
    %                        respect to the rotated projectors
    %            P0_/P1_: The rotated projectors

    % Pauli matrices
    paulis = {[0 1; 1 0], [0 -1i; 1i 0], [1 0; 0 -1]};

    % Normalize the rotation axis
    normal_vector_norm = NormalVector/(norm(NormalVector));

    % Construct the matrix 1/2 * n . sigma
    sin_part = zeros(2,2);
    for i = 1:3
        sin_part = sin_part + paulis{i}*0.5*normal_vector_norm(i);
    end

    % Two-qubit system, extend the rotation to the full system by acting
    % with the single-qubit rotation on the target qubit only
    if Qubits == 2
        if TargetQubit == 1
            rotation_matrix = kron(cos(Angles/2)*eye(2) - 2i*sin(Angles/2)*sin_part, eye(2));
        elseif TargetQubit == 2
            rotation_matrix = kron(eye(2), cos(Angles/2)*eye(2) - 2i*sin(Angles/2)*sin_part);
        end

    % Three-qubit system, extend the rotation to the full system by acting
    % with the single-qubit rotation on the target qubit only
    elseif Qubits == 3
        if TargetQubit == 1
            rotation_matrix = kron(cos(Angles/2)*eye(2) - 2i*sin(Angles/2)*sin_part, eye(4));
        elseif TargetQubit == 2
            rotation_matrix = kron(eye(2), kron(cos(Angles/2)*eye(2) - 2i*sin(Angles/2)*sin_part, eye(2)));
        elseif TargetQubit == 3
            rotation_matrix = kron(eye(4), cos(Angles/2)*eye(2) - 2i*sin(Angles/2)*sin_part);
        end
    end

    % Compute the rotated set of projectors
    P0_ = rotation_matrix*P0*rotation_matrix';
    P1_ = rotation_matrix*P1*rotation_matrix';

    % Extract the classical structure corresponding to the rotated
    % projectors
    RhoSigmaCl = P0_*RhoSigma*P0_ + P1_*RhoSigma*P1_;

    % Determine the overlap with the separable structure (shifted by the
    % identity contribution for stability), only consider real part
    overlap = real(trace(RhoSigmaCl*(RhoSigma + 0.25*eye(2^Qubits))));

    % Maximizing the overlap (OptStructure = 1) corresponds to minimizing
    % the cost, for OptStructure = -1 the overlap is minimized
    Cost = OptStructure*(1 - overlap);
end
