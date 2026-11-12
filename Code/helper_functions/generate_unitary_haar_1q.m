% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function U = generate_unitary_haar_1q()

    %!!! Generate a Haar-random single-qubit (2x2) unitary matrix, i.e., a
    %    unitary which is uniformly distributed with respect to the Haar
    %    measure on U(2). The unitary is obtained from the QR decomposition
    %    of a complex Ginibre matrix (Gaussian random matrix), where the
    %    phases of the diagonal of R are corrected to ensure the correct
    %    distribution (Mezzadri's algorithm)
    %
    %    input: none
    %
    %    output: U: The Haar-random 2x2 unitary matrix

    % Complex Gaussian random matrix (Ginibre ensemble)
    z = (randn(2) + 1i * randn(2)) / sqrt(2);

    % QR decomposition, Q is unitary but not yet Haar distributed
    [Q, R] = qr(z);

    % Normalize the phases of the diagonal elements of R and multiply them
    % column-wise onto Q to make the distribution Haar
    d = diag(R);
    ph = d ./ abs(d);
    U = Q .* ph.';
end
