% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [VNE] = get_von_neumann_entropy(DensityOp)

    %!!! Compute the von Neumann entropy of a given density operator
    %
    %    input: DensityOp: Density operator, must be Hermitian,
    %    positive-semidefinite
    %
    %    output: VNE: The von Neumann entropy
    
    % Check input parameters
    arguments
        DensityOp (:,:) {mustBeNumeric}
    end

    % Validate that Density operator is Hermitian and positive-semidefinite
    tol = 1e-8;
    isHermitian = norm(DensityOp - DensityOp', 'fro') <= tol * max(1, norm(DensityOp, 'fro'));

    % Is the input operator a valid density operator?
    if isHermitian
        lambdaMin = min(eig(DensityOp));
        isPSD = lambdaMin >= -tol * max(1, norm(DensityOp, 'fro'));
        if ~isPSD
            error("The input matrix is not a valid density operator. " + ...
                  "Ensure the matrix to be positive-semidefinite");
        end
    else
        error("The input matrix is not a valid density operator. " + ...
              "Ensure the matrix to be Hermitian");
    end

    % Compute the eigendecomposition of the density operator
    [~,D] = eig(DensityOp);
    eigvalues = diag(D);
    
    % Compute the von Neumann entropy
    VNE = 0;
    for ii = 1:size(eigvalues,1)
        % Continue if eigenvalue is smaller than numerical tolerance
        if abs(eigvalues(ii)) < tol
            continue;
        end

        VNE = VNE - eigvalues(ii) * log(eigvalues(ii));
    end
end


