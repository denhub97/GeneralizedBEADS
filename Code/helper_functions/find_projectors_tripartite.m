function [RhoSigmaClTri, RhoMClTri, RhoEClTri] = find_projectors_tripartite(DensityOpTri, RhoSigmaTri, RhoMTri, RhoETri, TargetQubits, OptStructure)

    %!!! Find the single-qubit measurement projectors on two target qubits
    %    of a 3-qubit state and apply them to the correlation structures of
    %    the state. The projectors are optimized for the density operator
    %    (maximal classical correlations, minimal conditional entropy of the
    %    remaining qubit). If the optimal orientations of the projectors
    %    are not unique, i.e., if the optimal measurement direction of a
    %    target qubit has a degree of freedom, the projectors are rotated
    %    such that the overlap with the separable structure of the state is
    %    maximal (OptStructure = 1) or minimal (OptStructure = -1), under
    %    the constraint that the optimal conditional entropy is maintained.
    %    The degree of freedom is determined from repeated optimizations
    %    (unique orientation, rotation about an axis, or arbitrary
    %    orientation) and whether the angles between the two projector
    %    directions are correlated. This analysis and the subsequent
    %    rotation of the projectors are only performed if OptStructure is
    %    not 0, otherwise the optimized projectors found first are applied
    %    to the structures directly
    %
    %    input: DensityOpTri: The density operator describing the 3-qubit
    %                         state (8x8)
    %           RhoSigmaTri: The separable structure of the density
    %                        operator
    %           RhoMTri: The mutual structure of the density operator
    %           RhoETri: The connected (entanglement) structure of the
    %                    density operator (an empty matrix can be passed,
    %                    it is then replaced by a zero matrix)
    %           TargetQubits: The two qubits which are measured
    %           OptStructure: Structure which the projectors are optimized
    %                            for (integer): -1 = entangling structure
    %                            (the overlap with the separable structure
    %                            is minimized), 0 = no structure
    %                            optimization, 1 = separable structure (the
    %                            overlap with the separable structure is
    %                            maximized). For -1 and 1 the degree of
    %                            freedom of the optimal projector
    %                            orientations is determined and used
    %
    %    output: RhoSigmaClTri: The classical structure of the separable
    %                           structure, i.e., RhoSigmaTri after the
    %                           measurement
    %            RhoMClTri: The mutual structure after the measurement
    %            RhoEClTri: The connected structure after the measurement

    % Handle default parameters
    if isempty(RhoETri)
        RhoETri = zeros(size(DensityOpTri,1));
    end

    % Define Paulis
    pauli_x = [0 1; 1 0];
    pauli_y = [0 -1i; 1i 0];
    pauli_z = [1 0; 0 -1];

    % Add identity to DensityOpTri to ensure positive semidefiniteness
    [~,D] = eig(DensityOpTri);
    min_eig = min(real(diag(D)));

    rho_plusid = (DensityOpTri+eye(8)*abs(min_eig))/trace(DensityOpTri+eye(8)*abs(min_eig));

    tolerance = 1e-4;
    runs = 10;
    bloch_vectors1 = zeros(runs,3);
    bloch_vectors2 = zeros(runs,3);
    angles_betweenprojectors = zeros(runs,1);
    i = 1;
    min_cost = 2;

    % Repeat the optimization of the projectors to collect the Bloch
    % vectors of all optimal projector directions
    while i <= runs

        % Optimize classical mutual information
        [cost,P00,P01,P10,P11,P0_1,P0_2] = opt_classcorrel3_haar_seq(rho_plusid, TargetQubits);

        % Collect all projectors by their corresponding Bloch vectors
        % Define the Bloch vectors based on the projector

        % A better optimum has been found, optimize again
        if i == 1 && (min_cost - cost) > tolerance
            min_cost = min(min_cost,cost);
            continue;
        end

        min_cost = min(min_cost,cost);

        % Bloch vector 1
        n_x_1 = trace(P0_1*pauli_x);
        n_y_1 = trace(P0_1*pauli_y);
        n_z_1 = trace(P0_1*pauli_z);

        bloch_vector1 = [n_x_1,n_y_1,n_z_1];
        bloch_vectors1(i,:) = bloch_vector1;

        % Bloch vector 2
        n_x_2 = trace(P0_2*pauli_x);
        n_y_2 = trace(P0_2*pauli_y);
        n_z_2 = trace(P0_2*pauli_z);

        bloch_vector2 = [n_x_2,n_y_2,n_z_2];
        bloch_vectors2(i,:) = bloch_vector2;

        % Angle between the two projector directions
        angle_proj = acos((bloch_vector1*bloch_vector2')/(norm(bloch_vector1)*norm(bloch_vector2)));
        angles_betweenprojectors(i) = angle_proj;

        % Sanity check, reset if we found a better solution
        if (min_cost - cost) > tolerance
            bloch_vectors1 = zeros(runs,3);
            bloch_vectors2 = zeros(runs,3);
            angles_betweenprojectors = zeros(runs,1);
            i = 1;
            continue;
        end
        i = i+1;
    end

    % Only if requested, analyze the degree of freedom of the optimal
    % projectors
    if OptStructure ~= 0
        % Now we try to find the best overlap of the optimized CMI with all
        % separable correlations
    
        % Check whether the projectors lie on two poles
        min_polecost1 = 1;
        min_polecost2 = 1;
        optset = zeros(3,2);
        % Rotational modes for projectors = elements 1 and 2, angle correlation = 3
        optimmode = zeros(3,1);
    
        % Fit a pair of antipodal poles to the Bloch vectors of each target
        % qubit (5 attempts each), the best fit is kept
        for ii = 1:5
            [~, polecost1, the1, phi1] = optPoles(bloch_vectors1);
            if polecost1 < min_polecost1
                min_polecost1 = polecost1;
                bestcoord1 = [the1; phi1];
            end
    
            [~, polecost2, the2, phi2] = optPoles(bloch_vectors2);
            if polecost2 < min_polecost2
                min_polecost2 = polecost2;
                bestcoord2 = [the2; phi2];
            end
        end
    
        polethreshold = 1e-4;
        conethreshold = runs*1e-5;
    
        % Analyze projector on target 1
        if min_polecost1 < polethreshold
            disp("Found unique projectors for target 1");
            optimmode(1) = 1;
            optset(:,1) = [min_polecost1; bestcoord1];
        else
            % Additional check for more complex geometries
            min_conecost1 = 1;
    
            % Fit a double cone with the tips on the axis to the Bloch vectors
            h1_ = 10;
            for ii = 1:runs
                [~, h1, conecost1, the1, phi1] = optTip2TipCone(bloch_vectors1);
                if conecost1 < min_conecost1
                    min_conecost1 = conecost1;
                    bestcoord1 = [the1; phi1];
                    h1_ = h1;
                end
            end
    
            if abs(h1_) < 1 && abs(h1_) >= 0 && min_conecost1 < conethreshold
                disp("Found rotational degree of freedom of projectors for target 1");
                optimmode(1) = 2;
                optset(:,1) = [min_conecost1; bestcoord1];
            else
                disp("Found abitrary or complex distribution of projectors for target 1");
                optimmode(1) = 3;
                optset(:,1) = [min_conecost1; rand(); rand()];
            end
        end
    
        % Analyze projector on target 2
        if min_polecost2 < polethreshold
            disp("Found unique projectors for target 2");
            optimmode(2) = 1;
            optset(:,2) = [min_polecost2; bestcoord2];
        else
            % Additional check for more complex geometries
            min_conecost2 = 1;
    
            % Fit a double cone with the tips on the axis to the Bloch vectors
            h2_ = 10;
            for ii = 1:runs
                [~, h2, conecost2, the2, phi2] = optTip2TipCone(bloch_vectors2);
                if conecost2 < min_conecost2
                    min_conecost2 = conecost2;
                    bestcoord2 = [the2; phi2];
                    h2_ = h2;
                end
            end
    
            if abs(h2_) < 1 && abs(h2_) >= 0 && min_conecost2 < conethreshold
                disp("Found rotational degree of freedom of projectors for target 2");
                optimmode(2) = 2;
                optset(:,2) = [min_conecost2; bestcoord2];
            else
                disp("Found abitrary or complex distribution of projectors for target 2");
                optimmode(2) = 3;
                optset(:,2) = [min_conecost2; rand(); rand()];
            end
        end
    
        % Figure out whether the angles between projectors are correlated
        ang_threshold = 1e-4;
        angles_corr = angles_betweenprojectors;
        angles_corr(abs(angles_corr - pi) < ang_threshold) = 0;
        angles_corr(abs(angles_corr) < ang_threshold) = 0;
        angles_correlated = max(abs(angles_corr - angles_corr(1))) < ang_threshold;
    
        if angles_correlated
            disp("The angles between projectors are correlated");
            optimmode(3) = 1;
        else
            disp("The angles between projectors are not correlated");
            optimmode(3) = 0;
        end
    end

    % Special case, projector directions fixed or no analysis of the degree
    % of freedom requested (OptStructure 0), no further steps required
    if  OptStructure == 0 || (optimmode(1) == 1 && optimmode(2) == 1)
        RhoMClTri = P00*RhoMTri*P00 + ...
                    P01*RhoMTri*P01 + ...
                    P10*RhoMTri*P10 + ...
                    P11*RhoMTri*P11;

        RhoSigmaClTri = P00*RhoSigmaTri*P00 + ...
                        P01*RhoSigmaTri*P01 + ...
                        P10*RhoSigmaTri*P10 + ...
                        P11*RhoSigmaTri*P11;

        RhoEClTri = P00*RhoETri*P00 + ...
                    P01*RhoETri*P01 + ...
                    P10*RhoETri*P10 + ...
                    P11*RhoETri*P11;
    else

        % Flatten the optimization parameters
        optparams = optset(:);

        % Adjust lower and upper bounds
        lb = optparams;
        ub = optparams;

        lb(1) = 0;
        lb(4) = 0;
        ub(1) = 2*pi;
        ub(4) = 2*pi;

        if optimmode(1) == 3
            lb(2:3) = [0, 0];
            ub(2:3) = [pi, 2*pi];
        end

        if optimmode(2) == 3
            lb(5:6) = [0, 0];
            ub(5:6) = [pi, 2*pi];
        end

        % Initialize arbitrary rotation angles
        optparams(1) = 2*pi*rand();
        optparams(4) = 2*pi*rand();

        % Get the optimum conditional information
        opt_condinfo = calc_conditionalinfo3(rho_plusid, P00, P01, P10, P11, TargetQubits);
        epsilon = 0.01 * opt_condinfo;

        % Optimization setup for each fmincon
        Options = optimset('HessUpdate','lbfgs','TolX',1e-8,...
            'TolFun',1e-8,'MaxIter', 5000, 'display','iter', 'UseParallel',true,...
            'maxfuneval',1e10, 'LargeScale', 'off', 'GradObj', 'off', ...
            'DerivativeCheck', 'off');

        % Global optimization problem, maximize the overlap with the
        % separable structure while the conditional entropy stays optimal
        problem = createOptimProblem('fmincon', ...
        'objective', @(params) cost_sigma_cl3(params, RhoSigmaTri, P00, P01, P10, P11, optimmode(3), TargetQubits, OptStructure), ...
        'x0', optparams, ...
        'lb', lb, 'ub', ub, ...
        'nonlcon', @(params) check_cmi(params, rho_plusid, opt_condinfo, epsilon, P00, P01, P10, P11, optimmode(3), TargetQubits), ...
        'options', Options);

        gs = GlobalSearch('NumTrialPoints', 1000, 'Display', 'iter', 'NumStageOnePoints', 200);
        params_best = run(gs, problem);

        % Apply the optimally rotated projectors to the structures
        [~, RhoSigmaClTri, P00_, P01_, P10_, P11_] = cost_sigma_cl3(params_best, RhoSigmaTri, P00, P01, P10, P11, optimmode(3), TargetQubits, OptStructure);
        RhoMClTri = P00_*RhoMTri*P00_ + ...
                    P01_*RhoMTri*P01_ + ...
                    P10_*RhoMTri*P10_ + ...
                    P11_*RhoMTri*P11_;

        RhoEClTri = P00_*RhoETri*P00_ + ...
                    P01_*RhoETri*P01_ + ...
                    P10_*RhoETri*P10_ + ...
                    P11_*RhoETri*P11_;
    end
end

function [c, ceq] = check_cmi(params, rho, CMI_target, epsilon, P00, P01, P10, P11, ang_corr, targetqubits)

    % Nonlinear constraints for the optimization of the projector
    % rotation: the conditional entropy of the rotated projectors must stay
    % within the interval [CMI_target, CMI_target + epsilon]
    %
    %    input: params: The rotation parameters of the projectors
    %           rho: The density operator for which the conditional
    %                entropy is evaluated
    %           CMI_target: The optimal conditional entropy
    %           epsilon: The allowed deviation from the optimal value
    %           P00/P01/P10/P11: The initial projectors
    %           ang_corr: Are the projector orientations correlated?
    %           targetqubits: The two target qubits
    %
    %    output: c: The inequality constraints (c <= 0)
    %            ceq: The equality constraints (none)

    % Pauli matrices
    sigma_x = [0, 1; 1, 0];
    sigma_y = [0, -1i; 1i, 0];
    sigma_z = [1, 0; 0, -1];

    % Compute components of the rotation axis in Cartesian coordinates
    nx1 = sin(params(2)) * cos(params(3));
    ny1 = sin(params(2)) * sin(params(3));
    nz1 = cos(params(2));

    % Construct n . sigma
    n_dot_sigma1 = nx1 * sigma_x + ny1 * sigma_y + nz1 * sigma_z;

    % Compute the unitary rotation matrix
    U1 = cos(params(1)/2) * eye(2) - 1i * sin(params(1)/2) * n_dot_sigma1;

    % If the angles between the projection directions are not correlated
    if ang_corr == 0
        % Compute components of the rotation axis in Cartesian coordinates
        nx2 = sin(params(5)) * cos(params(6));
        ny2 = sin(params(5)) * sin(params(6));
        nz2 = cos(params(5));

        % Construct n dot sigma
        n_dot_sigma2 = nx2 * sigma_x + ny2 * sigma_y + nz2 * sigma_z;

        % Compute the unitary rotation matrix
        U2 = cos(params(4)/2) * eye(2) - 1i * sin(params(4)/2) * n_dot_sigma2;
    else
        U2 = U1;
    end

    % Sort the target qubit list and construct the rotation of the full
    % system, acting on the two target qubits only
    targetqubits_sorted = sort(targetqubits, 'ascend');
    if isequal(targetqubits_sorted, [1,2])
        U = kron(kron(U1,U2), eye(2));
    elseif isequal(targetqubits_sorted, [1,3])
        U = kron(kron(U1,eye(2)), U2);
    else
        U = kron(kron(eye(2),U1), U2);
    end

    % Compute the rotated set of projectors
    P00_ = U*P00*U';
    P01_ = U*P01*U';
    P10_ = U*P10*U';
    P11_ = U*P11*U';

    % Conditional entropy of the rotated projectors
    CMI_val = calc_conditionalinfo3(rho, P00_, P01_, P10_, P11_, targetqubits);

    % Inequality constraints
    c = [
         CMI_target - CMI_val;
         CMI_val - (CMI_target + epsilon);
        ];

    % No equality constraint
    ceq = [];
end
