% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [Cost, P0, P1] = opt_classcorrel3_haar_single(DensityOp, TargetQubit)

    %!!! Optimization function for finding the projectors of the measurement
    %    basis in which, by measuring a single target qubit of a tripartite
    %    state, the reduced density operator of the remaining qubits
    %    possesses the most information (minimally mixed = minimal
    %    entropy). The projector orientation is optimized starting from a
    %    Haar-random initial orientation. The conditional entropy computed
    %    by cost_classcorr3_haar_seq is minimized -> maximized classical
    %    correlation
    %
    %    input: DensityOp: The density operator describing the tripartite
    %                      state
    %           TargetQubit: The qubit which is measured
    %
    %    output: ClassCorr: The classical correlation cost value
    %                       (conditional entropy) of the optimized
    %                       measurement
    %            P0/P1: The optimized projectors on the measured qubit,
    %                   extended to the full 3-qubit system

    % Initialize 1 qubit unitary matrix vector for optimization
    % Sample a Haar-random unitary
    U1 = generate_unitary_haar_1q();
    [theta1, phi1, lambda1] = extract_euler_angles_1q(U1);

    % Corresponding Euler angles as initial optimization parameters
    U_ini = [theta1, phi1, lambda1];

    % Settings of the quasi-Newton optimizer
    Options = optimset('HessUpdate','lbfgs','TolX',1e-8,...
        'TolFun',1e-8,'MaxIter', 5000, 'display','iter', 'UseParallel',true,...
        'maxfuneval',1e10, 'LargeScale', 'off', 'GradObj', 'off', ...
        'DerivativeCheck', 'off');

    % Optimize the Euler angles, i.e., the projector orientation
    calc = @(U_coeff) cost_classcorr3_haar_seq(U_coeff,DensityOp, TargetQubit);
    [U_opt, Cost, exitflag, output] = fminunc(calc, U_ini, Options);

    % Evaluate the full-system projectors for the optimized coefficients
    [Cost, P0, P1, ~] = cost_classcorr3_haar_seq(U_opt, DensityOp, TargetQubit);
end
