% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [Theta, Phi, Lambda] = extract_euler_angles_1q(U)

    %!!! Extract the Euler angles [theta, phi, lambda] of a single-qubit
    %    (2x2) unitary matrix in the ZYZ convention, i.e., up to a global
    %    phase,
    %
    %        U = [cos(theta/2),            -e^{i*lambda}*sin(theta/2);
    %             e^{i*phi}*sin(theta/2),   e^{i*(phi+lambda)}*cos(theta/2)]
    %
    %    which corresponds to Rz(phi) * Ry(theta) * Rz(lambda). This is
    %    the parametrization used, e.g., in cost_classcorr3_haar_seq. The
    %    extraction is independent of the global phase of the unitary,
    %    since only phase differences of the matrix entries are used
    %
    %    input: U: The 2x2 unitary matrix
    %
    %    output: Theta: The polar Euler angle, in [0, pi]
    %            Phi: The first azimuthal Euler angle, in [0, 2*pi)
    %            Lambda: The second azimuthal Euler angle, in [0, 2*pi)

    % Ensure U is 2x2 unitary
    if ~isequal(size(U), [2, 2])
        error('Input must be a 2x2 unitary matrix.');
    end

    % Calculate theta from the magnitudes of the first column, |U(1,1)| =
    % cos(theta/2) and |U(2,1)| = sin(theta/2). The atan2 formulation is
    % numerically robust, in contrast to acos close to 1, and yields
    % theta in [0, pi]
    Theta = 2 * atan2(abs(U(2,1)), abs(U(1,1)));

    % Diagonal unitary, only the sum phi + lambda is defined
    if abs(U(2,1)) < 1e-12
        Phi = 0;
        Lambda = angle(U(2,2)) - angle(U(1,1));  % Arbitrary
    % Anti-diagonal unitary, only the difference phi - lambda is defined
    elseif abs(U(1,1)) < 1e-12
        Lambda = 0;
        Phi = angle(U(2,1)) - angle(-U(1,2));  % Arbitrary
    % General case, phases relative to the phase of the first entry
    else
        Phi = angle(U(2,1)) - angle(U(1,1));
        Lambda = angle(-U(1,2)) - angle(U(1,1));
    end

    % Wrap angles into [0, 2*pi)
    Theta = mod(Theta, 2*pi);
    Phi = mod(Phi, 2*pi);
    Lambda = mod(Lambda, 2*pi);
end
