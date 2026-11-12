% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [ClassCorr, P00, P01, P10, P11, P0_1_, P0_2_] = opt_classcorrel3_haar_seq(DensityOpTri, TargetQubits)

    %!!! Optimization function for finding the projectors of the measurement
    %    basis in which, by measuring two target qubits of a tripartite
    %    state, the reduced density operator of the remaining qubit
    %    possesses the most information (minimally mixed = minimal
    %    entropy). The single-qubit projectors on the two target qubits are
    %    optimized sequentially, starting from Haar-random initial
    %    orientations. The conditional entropy computed by
    %    cost_classcorr3_haar_seq is minimized -> maximized classical
    %    correlation
    %
    %    input: DensityOpTri: The density operator describing the tripartite
    %                         state
    %           TargetQubits: The list of the two target qubits which are
    %                         measured
    %
    %    output: ClassCorr: The classical correlation cost value (joint
    %                       conditional entropy) of the optimized
    %                       measurements
    %            P00/P01/P10/P11: Full set of joint projectors on the two
    %                             target qubits, extended to the full system
    %            P0_1_/P0_2_: The single-qubit projectors corresponding to
    %                         outcome 0 on the first and second target qubit

    % Sort the list of target qubits in ascending order, the first target
    % qubit is optimized first
    TargetQubits = sort(TargetQubits,"ascend");

    % Initialize 1 qubit unitary matrix vector for optimization
    % Sample two Haar-random unitaries
    U1 = generate_unitary_haar_1q();
    [theta1, phi1, lambda1] = extract_euler_angles_1q(U1);

    U2 = generate_unitary_haar_1q();
    [theta2, phi2, lambda2] = extract_euler_angles_1q(U2);

    % Corresponding Euler angles as initial optimization parameters
    U_ini = [theta1, phi1, lambda1];
    U_ini2 = [theta2, phi2, lambda2];

    % Settings of the quasi-Newton optimizer
    Options = optimset('HessUpdate','lbfgs','TolX',1e-8,...
        'TolFun',1e-8,'MaxIter', 5000, 'display','iter', 'UseParallel',true,...
        'maxfuneval',1e10, 'LargeScale', 'off', 'GradObj', 'off', ...
        'DerivativeCheck', 'off');

    % Optimize the projector orientation of the first target qubit
    calc = @(U_coeff) cost_classcorr3_haar_seq(U_coeff,DensityOpTri,TargetQubits(1));
    [U_opt, cost_final, exitflag, output] = fminunc(calc, U_ini, Options);

    % Evaluate the projectors of the first target qubit for the optimized
    % coefficients
    [ClassCorr,P0_1,P1_1, P0_1_] = cost_classcorr3_haar_seq(U_opt,DensityOpTri,TargetQubits(1));

    % Post-measurement state after measuring the first target qubit
    % (outcomes not resolved)
    rho_p = P0_1*DensityOpTri*P0_1 + ...
            P1_1*DensityOpTri*P1_1;

    % Optimize the projector orientation of the second target qubit for the
    % post-measurement state
    calc2 = @(U_coeff2) cost_classcorr3_haar_seq(U_coeff2,rho_p,TargetQubits(2));
    [U_opt2, cost_final2, exitflag, output] = fminunc(calc2, U_ini2, Options);

    % Evaluate the projectors of the second target qubit for the optimized
    % coefficients
    [ClassCorr,P0_2,P1_2, P0_2_] = cost_classcorr3_haar_seq(U_opt2,rho_p,TargetQubits(2));

    % Combine to the full set of joint projectors on both target qubits
    P00 = P0_1*P0_2;
    P01 = P0_1*P1_2;
    P10 = P1_1*P0_2;
    P11 = P1_1*P1_2;

    % Sort the list of target qubits in descending order for efficient
    % partial trace operations
    TargetQubits = sort(TargetQubits,"descend");

    % Collect projectors, probabilities and outcomes
    projectors = {P00, P01, P10, P11};
    probabilities = zeros(1,4);
    outcomes = cell(1,4);

    % Compute the probability and the post-measurement state of the
    % remaining qubit for each joint outcome
    for pp = 1:numel(probabilities)
        % The current projector
        P = projectors{pp};

        probabilities(pp) = trace(P*DensityOpTri);
        outcomes{pp} = partial_trace(partial_trace(P*DensityOpTri*P,TargetQubits(1)),TargetQubits(2)) / probabilities(pp);
    end

    % Calcuate the conditional entropy in terms of the outcome
    % eigenvalue decomposition, eigenvalues smaller than tolerance do not
    % contribute to the entropy (evaluate to 0). This circumvents the
    % problem of logm() = NaN, when rho does not have full rank.

    % The numerical tolerance for eigenvalues
    tolerance = 1e-8;

    ClassCorr = 0;
    for ii = 1:numel(probabilities)
        p = probabilities(ii);
        outcome = outcomes{ii};

        % Compute the eigendecomposition of the outcome
        [~,D] = eig(outcome);
        eigvalues = diag(D);

        % Compute the von Neumann entropy of the outcome
        S = 0;
        for jj = 1:size(eigvalues,1)
            % Continue if eigenvalue is smaller than tolerance
            if abs(eigvalues(jj)) < tolerance
                continue;
            end

            S = S + eigvalues(jj)*log(eigvalues(jj));
        end

        % Weight with the outcome probability
        S = p*(-S);

        % Add to the joint conditional entropy, only consider real part
        ClassCorr = ClassCorr + real(S);
    end
end
