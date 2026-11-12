function [RedDensityOp] = partial_trace(DensityOp, Subsystem)

    %!!! Compute the partial trace of a density operator of a composite
    %    system of up to 4 qubits. The partial trace is taken over the
    %    specified qubit by summing over a basis of this qubit, resulting
    %    in the reduced density operator of the remaining qubits
    %
    %    input: DensityOp: The density operator of the composite system
    %                      (4x4, 8x8, or 16x16 for 2, 3, or 4 qubits)
    %           Subsystem: The qubit which the partial trace is taken over
    %                      (1 to the number of qubits)
    %
    %    output: RedDensityOp: The reduced density operator of the
    %                          remaining qubits

    % Check input parameters
    arguments
        DensityOp (:,:) {mustBeNumeric}
        Subsystem (1,1) {mustBeInteger, mustBePositive}
    end

    % The density operator must be square and describe 2, 3, or 4 qubits
    assert(size(DensityOp,1) == size(DensityOp,2) && ismember(size(DensityOp,1), [4 8 16]), ...
           "partial_trace:invalidDimension", ...
           "DensityOp must be a square matrix of size 4, 8, or 16, but is %d x %d.", ...
           size(DensityOp,1), size(DensityOp,2));

    % The subsystem must be one of the qubits of the system
    assert(Subsystem <= log2(size(DensityOp,1)), ...
           "partial_trace:invalidSubsystem", ...
           "Subsystem must be between 1 and %d, but is %d.", ...
           log2(size(DensityOp,1)), Subsystem);

    % Computational basis of a single qubit
    basis_vecs = {[1;0],[0;1]};

    % Two-qubit system
    if size(DensityOp,1) == 4

        RedDensityOp = zeros(2);

        % Trace out the first qubit
        if Subsystem == 1
            for j = 1:size(basis_vecs,2)
                RedDensityOp = RedDensityOp+kron(basis_vecs{j}',eye(2))*DensityOp*kron(basis_vecs{j},eye(2));
            end

        % Trace out the second qubit
        elseif Subsystem == 2
            for j = 1:size(basis_vecs,2)
                RedDensityOp = RedDensityOp+kron(eye(2),basis_vecs{j}')*DensityOp*kron(eye(2),basis_vecs{j});
            end
        end

    % Three-qubit system
    elseif size(DensityOp,1) == 8

        RedDensityOp = zeros(4);

        % Trace out the first qubit
        if Subsystem == 1
            for j = 1:size(basis_vecs,2)
                RedDensityOp = RedDensityOp+kron(basis_vecs{j}',eye(4))*DensityOp*kron(basis_vecs{j},eye(4));
            end

        % Trace out the second qubit
        elseif Subsystem == 2
            for j = 1:size(basis_vecs,2)
                RedDensityOp = RedDensityOp+kron(kron(eye(2),basis_vecs{j}'),eye(2))*DensityOp*kron(kron(eye(2),basis_vecs{j}),eye(2));
            end

        % Trace out the third qubit
        elseif Subsystem == 3
            for j = 1:size(basis_vecs,2)
                RedDensityOp = RedDensityOp+kron(eye(4),basis_vecs{j}')*DensityOp*kron(eye(4),basis_vecs{j});
            end
        end

    % Four-qubit system
    elseif size(DensityOp,1) == 16

        RedDensityOp = zeros(8);

        % Trace out the first qubit
        if Subsystem == 1
            for j = 1:size(basis_vecs,2)
                RedDensityOp = RedDensityOp+kron(basis_vecs{j}',eye(8))*DensityOp*kron(basis_vecs{j},eye(8));
            end

        % Trace out the second qubit
        elseif Subsystem == 2
            for j = 1:size(basis_vecs,2)
                RedDensityOp = RedDensityOp+kron(kron(eye(2),basis_vecs{j}'),eye(4))*DensityOp*kron(kron(eye(2),basis_vecs{j}),eye(4));
            end

        % Trace out the third qubit
        elseif Subsystem == 3
            for j = 1:size(basis_vecs,2)
                RedDensityOp = RedDensityOp+kron(kron(eye(4),basis_vecs{j}'),eye(2))*DensityOp*kron(kron(eye(4),basis_vecs{j}),eye(2));
            end

        % Trace out the fourth qubit
        elseif Subsystem == 4
            for j = 1:size(basis_vecs,2)
                RedDensityOp = RedDensityOp+kron(eye(8),basis_vecs{j}')*DensityOp*kron(eye(8),basis_vecs{j});
            end
        end
    end
end
