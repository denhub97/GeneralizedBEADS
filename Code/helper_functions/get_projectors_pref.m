function [P0, P1] = get_projectors_pref(RAxis, TargetQubit, N)

    %!!! Compute the projectors of a single-qubit measurement along a given
    %    (preferred) axis, extended to the full 2-qubit or 3-qubit system.
    %    The projectors are P0 = (1 + n . sigma)/2 and P1 = (1 - n . sigma)/2
    %    for the unit vector n = RAxis, i.e., the projectors onto the
    %    states with spin up and spin down along the measurement axis, and
    %    act on the target qubit only
    %
    %    input: RAxis: The measurement direction as unit vector
    %                  [x; y; z] on the Bloch sphere
    %           TargetQubit: The qubit which is measured
    %           N: The number of qubits of the system (2 or 3, default 2)
    %
    %    output: P0/P1: The projectors for the two measurement outcomes,
    %                   extended to the full system

    % Handle default parameters
    if nargin < 3
        N = 2;
    end

    % Check that the measurement axis is a three-dimensional unit vector
    assert(isnumeric(RAxis) && numel(RAxis) == 3 && abs(norm(RAxis) - 1) < 1e-10, ...
           "get_projectors_pref:invalidAxis", ...
           "RAxis must be a numeric unit vector with 3 components.");

    % Check that the system size is 2 or 3 qubits and that the target
    % qubit is one of its qubits
    assert(isscalar(N) && ismember(N, [2 3]), ...
           "get_projectors_pref:invalidNumberOfQubits", ...
           "N must be 2 or 3.");
    assert(isscalar(TargetQubit) && ismember(TargetQubit, 1:N), ...
           "get_projectors_pref:invalidTargetQubit", ...
           "TargetQubit must be an integer between 1 and N = %d.", N);

    % Single-qubit projector onto the measurement direction
    P = 1/2 * [1+RAxis(3), RAxis(1)-1i*RAxis(2);
               RAxis(1)+1i*RAxis(2), 1-RAxis(3)];

    % Single-qubit projector onto the opposite direction
    P_ = 1/2 * [1-RAxis(3), -RAxis(1)+1i*RAxis(2);
                -RAxis(1)-1i*RAxis(2), 1+RAxis(3)];

    % Two-qubit system, extend the projectors to the full system by acting
    % on the target qubit only
    if N == 2
        if TargetQubit == 1
            P0 = kron(P,eye(2));
            P1 = kron(P_,eye(2));
        elseif TargetQubit == 2
            P0 = kron(eye(2),P);
            P1 = kron(eye(2),P_);
        end

    % Three-qubit system, extend the projectors to the full system by
    % acting on the target qubit only
    elseif N == 3
        if TargetQubit == 1
            P0 = kron(kron(P, eye(2)), eye(2));
            P1 = kron(kron(P_, eye(2)), eye(2));
        elseif TargetQubit == 2
            P0 = kron(kron(eye(2), P), eye(2));
            P1 = kron(kron(eye(2), P_), eye(2));
        elseif TargetQubit == 3
            P0 = kron(kron(eye(2), eye(2)), P);
            P1 = kron(kron(eye(2), eye(2)), P_);
        end
    end
end
