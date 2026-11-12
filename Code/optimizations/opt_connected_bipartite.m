% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [EntNorm, States] = opt_connected_bipartite(DensityOp, Padd, PaddVectors)

    %!!! Optimization function for finding the pure state decomposition of
    %    a density operator with the smallest connected correlation
    %    function norm. The coefficients of a unitary mixing the
    %    eigenstates of the density operator are optimized, using the
    %    Frobenius norm of the connected correlations as cost function (see
    %    connected_norm_bipartite). Optionally, the decomposition can be
    %    padded with additional zero states
    %
    %    input: DensityOp: The density operator describing the state
    %           Padd: Shall the decomposition be padded with additional
    %                 (zero) states (boolean, default false)
    %           PaddVectors: Number of additional states for padding
    %                        (default 0)
    %
    %    output: EntNorm: The connected correlation function norm of the
    %                     optimized decomposition
    %            States: The cell array of pure states (including their
    %                    weights) constituting the optimal decomposition

    % Handle default parameters
    if nargin < 2
        Padd = false;
        PaddVectors = 0;
    end

    dim = size(DensityOp,1);

    % Start with an eigendecomposition of rho in order to have a convenient
    % starting point (basis vectors)
    [eigenv,eigenvals] = eig(DensityOp);
    eigenvalues = diag(eigenvals);

    % Find non-zero eigenvalues and their indices in order to have the
    % number of basis states equal the rank of rho
    tolerance = 10e-8;
    non_zero_indices = find(eigenvalues>tolerance);

    % Concatenate eigenvectors corresponding to non-zero eigenvalues
    nonzero_eigenvectors = [eigenv(:,non_zero_indices)];
    nonzero_eigenvalues = [eigenvalues(non_zero_indices)];
    rank_rho = length(nonzero_eigenvalues);

    % Create eigenstates by multiplying the square root of the eigenvalue
    % with the corresponding eigenvector
    eigenstates = nonzero_eigenvectors .* sqrt(nonzero_eigenvalues)';

    % Pad the set of eigenstates with zero vectors to allow decompositions
    % with more pure states than the rank of the density operator
    if Padd
        for ii = rank_rho+1 : rank_rho+PaddVectors
            if PaddVectors <= 0
                break;
            end
            eigenstates(:,ii) = 0;
        end
        rank_rho = rank_rho+PaddVectors;
    end

    % Randomly initialize the coefficients for the unitary
    U_ini = rand(1,rank_rho^2)';

    % Settings of the quasi-Newton optimizer (high accuracy)
    Options = optimset('HessUpdate','lbfgs','TolX',1e-12,...
        'TolFun',1e-12,'MaxIter', 5000, 'display','iter', 'UseParallel',true,...
        'MaxFunEvals',1e10, 'LargeScale', 'off', 'GradObj', 'off', ...
        'DerivativeCheck', 'off');

    % Optimize the unitary coefficients, i.e., the decomposition
    calc = @(U_coeff) connected_norm_bipartite(U_coeff,eigenstates,rank_rho,dim);
    [U_opt, cost_final, exitflag, output] = fminunc(calc, U_ini, Options);

    % Evaluate the connected correlation norm and the states for the
    % optimized coefficients
    [EntNorm,States] = connected_norm_bipartite(U_opt,eigenstates,rank_rho,dim);

    % Set small values to zero and normalize phase wrt biggest value
    tolerance = 1e-8;

    for k = 1:size(States,1)
        % Set entries smaller than tolerance to zero
        solutionk = States{k};
        solutionk(abs(solutionk)<tolerance) = 0;

        % Find the entry with the largest absolute value
        [max_value, idx] = max(abs(solutionk));
        max_element = solutionk(idx);

        % Calculate the phase of the maximum element
        max_element_phase = angle(max_element);

        % Adjust the phase of the entire vector to set the phase of max element to zero
        solutionk = solutionk * exp(-1i * max_element_phase);

        States{k} = solutionk;
    end

    % Sort the position of the arrays within the cell array States for
    % better compatibility with other final states from other optimizations

    % Create a matrix to store the real parts of the first entries and their indices
    first_entry_real_parts = zeros(size(States, 1), 2);

    % Extract the real parts of the first entries
    for k = 1:size(States, 1)
        solutionk = States{k};
        first_entry_real_parts(k, :) = [real(solutionk(1)), k];
    end

    % Sort based on the real parts in descending order
    sorted_entries = sortrows(first_entry_real_parts, -1);

    % Create a new cell array to store the sorted States
    sorted_final_states = cell(size(States));

    % Reorder the States based on the sorted indices
    for k = 1:size(States, 1)
        sorted_final_states{k} = States{sorted_entries(k, 2)};
    end

    % Update the original States with sorted States
    States = sorted_final_states;
end
