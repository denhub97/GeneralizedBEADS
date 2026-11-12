% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [Beads, EntNorm] = GetBeads3Qubits(PureStates123, PureStates12, PureStates13, PureStates23, BeadTypes, TargetQubits, GridPoints, ClassCorrRecursive, OptStructure, MeasurementDirection)

    %!!! Compute the spherical function values of all beads of a 3-qubit
    %    state, given as pure state decompositions of the density operator
    %    and of its three two-qubit marginals. The beads of the pure states
    %    are weighted with their squared norm and summed. Next to the beads
    %    of the total correlation (T), the connected correlations (E) and
    %    the compound part (C), the beads of the separable structure + its
    %    classical and nonclassical part, of the mutual structure and its
    %    classical part, and of the classical part of the connected
    %    structure are computed (using the measurement projectors found
    %    for the target qubits). The beads are computed for the full
    %    3-qubit state (single-qubit and trilinear beads) and for the
    %    bilinear beads of the marginals of qubits 1-2, 1-3, and 2-3. The
    %    tripartite entanglement norm, the two-tangles and the three-tangle
    %    of the state are computed as well
    %
    %    input: PureStates123: Cell array of the pure states constituting
    %                          the decomposition of the 3-qubit state
    %                          (the weights are contained in the state
    %                          vectors, i.e., sqrt(p)*|state>)
    %           PureStates12/PureStates13/PureStates23: Cell arrays of the
    %                          pure states constituting the decompositions
    %                          of the two-qubit marginals
    %           BeadTypes: List of the bead types (strings or character
    %                      arrays, at least one element)
    %           TargetQubits: The target qubits for the classical
    %                         structures, [q12, q13, q23, q1, q2] with the
    %                         target qubit of the marginals 1-2, 1-3, and
    %                         2-3 and the two target qubits of the
    %                         tripartite structure
    %           GridPoints: The beads resolution [theta, phi]
    %           ClassCorrRecursive: Shall the projectors of the classical
    %                               structures be determined recursively
    %                               for the marginals and the tripartite
    %                               structure separately (boolean,
    %                               default true)
    %           OptStructure: Structure which the classical beads are
    %                            optimized for (integer): -1 = entangling
    %                            structure (the overlap with the separable
    %                            structure is minimized), 0 = no structure
    %                            optimization, 1 = separable structure (the
    %                            overlap with the separable structure is
    %                            maximized) (default 0). It is passed to
    %                            find_projectors_tripartite and only used if
    %                            ClassCorrRecursive is true
    %           MeasurementDirection: Measurement direction [x y z] used
    %                            by find_projectors_bipartite for all
    %                            marginals if the optimal orientation is
    %                            arbitrary. If empty (default), the
    %                            direction is requested in the command window
    %
    %    output: Beads: Cell array with the beads, in the order: identity
    %                   bead, Q-beads 1 to 3, the even beads (bilinear
    %                   beads of the marginals 1-2, 1-3, 2-3, and the
    %                   trilinear beads tau2, tau3, tau4; each with T, E,
    %                   C, the separable structure, its classical and
    %                   nonclassical part, the mutual structure, and the
    %                   classical mutual structure), the odd beads (same
    %                   order, with the trilinear beads tau2, tau3, tau1),
    %                   the classical slice of connected beads, the target
    %                   qubits, the two-tangles, and the three-tangle
    %            EntNorm: The tripartite entanglement norm of the state

    % Handle default parameters
    if nargin < 8
        ClassCorrRecursive = true;
        OptStructure = 0;
        MeasurementDirection = [];
    elseif nargin < 9
        OptStructure = 0;
        MeasurementDirection = [];
    elseif nargin < 10
        MeasurementDirection = [];
    end

    % Check that the bead types are a non-empty list of strings or
    % character arrays
    assert((isstring(BeadTypes) || ischar(BeadTypes) || iscellstr(BeadTypes)) && ~isempty(BeadTypes), ...
           "GetBeads3Qubits:invalidBeadTypes", ...
           "BeadTypes must be a non-empty list of strings or character arrays.");

    % Determine which beads are required for the requested bead types
    BeadTypes = string(BeadTypes);

    % The separable structure (Sigma beads) is required for the mutual, the
    % separable, and all (non)classical beads
    NeedSigma = any(ismember(BeadTypes, ["M", "Sigma", "M_Cl", "Sigma_Cl", "D", "Sigma_D"]));

    % The measurement projectors and the structures after the measurement
    % are required for the classical beads and the nonclassical beads
    NeedProj = any(ismember(BeadTypes, ["M_Cl", "Sigma_Cl", "D", "Sigma_D"]));

    % The mutual structure (M beads) is required for the mutual beads and
    % for the projector analysis
    NeedM = any(ismember(BeadTypes, ["M", "M_Cl", "Sigma_Cl", "D", "Sigma_D"]));

    % Grid points for the beads
    res_theta = GridPoints(1);
    res_phi = GridPoints(2);

    % Lisa basis decomposition and mapping for 2 and 3 qubits
    % Load definitions of real Lisa tensors and (real and complex) spherical harmonics
    SH_Defs_Real;
    struc1 = load('LiSAReal.mat');
    LiSAReal = struc1.LiSAReal;
    struc2 = load('LiSAScalings.mat');
    LiSAScalings = struc2.LiSAScalings;

    % How many pure basis states contribute to mixed quantum state
    states_123 = size(PureStates123,1);
    dim_123 = size(PureStates123{1},1);

    states_12 = size(PureStates12,1);
    dim_12 = size(PureStates12{1},1);

    states_13 = size(PureStates13,1);
    dim_13 = size(PureStates13{1},1);

    states_23 = size(PureStates23,1);
    dim_23 = size(PureStates23{1},1);

    % Initialize beads_lists for collecting contributions from all pure states
    % and weigh them by norm squared of pure states

    % All beads, total correlation
    id_bead_mixed = zeros(size(Y00));
    bead1_mixed = zeros(size(Y00));
    bead2_mixed = zeros(size(Y00));
    bead3_mixed = zeros(size(Y00));

    bilin12_even_mixed = zeros(size(Y00));
    bilin12_odd_mixed = zeros(size(Y00));
    bilin13_even_mixed = zeros(size(Y00));
    bilin13_odd_mixed = zeros(size(Y00));
    bilin23_even_mixed = zeros(size(Y00));
    bilin23_odd_mixed = zeros(size(Y00));

    trilin_tau1_odd_mixed = zeros(size(Y00));
    trilin_tau2_even_mixed = zeros(size(Y00));
    trilin_tau2_odd_mixed = zeros(size(Y00));
    trilin_tau3_even_mixed = zeros(size(Y00));
    trilin_tau3_odd_mixed = zeros(size(Y00));
    trilin_tau4_mixed = zeros(size(Y00));

    % Bilinear terms, connected

    bilin12_even_mixed_E = zeros(size(Y00));
    bilin12_odd_mixed_E = zeros(size(Y00));
    bilin13_even_mixed_E = zeros(size(Y00));
    bilin13_odd_mixed_E = zeros(size(Y00));
    bilin23_even_mixed_E = zeros(size(Y00));
    bilin23_odd_mixed_E = zeros(size(Y00));

    trilin_tau1_odd_mixed_E = zeros(size(Y00));
    trilin_tau2_even_mixed_E = zeros(size(Y00));
    trilin_tau2_odd_mixed_E = zeros(size(Y00));
    trilin_tau3_even_mixed_E = zeros(size(Y00));
    trilin_tau3_odd_mixed_E = zeros(size(Y00));
    trilin_tau4_mixed_E = zeros(size(Y00));

    % Bilinear terms, compound
    bilin12_even_mixed_C = zeros(size(Y00));
    bilin12_odd_mixed_C = zeros(size(Y00));
    bilin13_even_mixed_C = zeros(size(Y00));
    bilin13_odd_mixed_C = zeros(size(Y00));
    bilin23_even_mixed_C = zeros(size(Y00));
    bilin23_odd_mixed_C = zeros(size(Y00));

    trilin_tau1_odd_mixed_C = zeros(size(Y00));
    trilin_tau2_even_mixed_C = zeros(size(Y00));
    trilin_tau2_odd_mixed_C = zeros(size(Y00));
    trilin_tau3_even_mixed_C = zeros(size(Y00));
    trilin_tau3_odd_mixed_C = zeros(size(Y00));
    trilin_tau4_mixed_C = zeros(size(Y00));

    rho_compound_mixed_123 = zeros(dim_123, dim_123);
    rho_connected_mixed_123 = zeros(dim_123, dim_123);

    rho_compound_mixed_12 = zeros(dim_12, dim_12);
    rho_connected_mixed_12 = zeros(dim_12, dim_12);

    rho_compound_mixed_13 = zeros(dim_13, dim_13);
    rho_connected_mixed_13 = zeros(dim_13, dim_13);

    rho_compound_mixed_23 = zeros(dim_23, dim_23);
    rho_connected_mixed_23 = zeros(dim_23, dim_23);
    EntNorm = 0;

    % Put together original rho for further analysis
    rho_original = zeros(dim_123,dim_123);
    rho_12 = zeros(dim_12, dim_12);
    rho_13 = zeros(dim_13, dim_13);
    rho_23 = zeros(dim_23, dim_23);

    % Beads of the separable structure, classical and nonclassical beads, mutual
    % beads, classical mutual and connected beads, these are only computed if
    % required and are initialized with zeros
    bilin12_even_Sigma = zeros(size(Y00));
    bilin12_even_SigmaCL = zeros(size(Y00));
    bilin12_even_SigmaD = zeros(size(Y00));
    bilin12_even_M = zeros(size(Y00));
    bilin12_even_MCl = zeros(size(Y00));
    bilin13_even_Sigma = zeros(size(Y00));
    bilin13_even_SigmaCL = zeros(size(Y00));
    bilin13_even_SigmaD = zeros(size(Y00));
    bilin13_even_M = zeros(size(Y00));
    bilin13_even_MCl = zeros(size(Y00));
    bilin23_even_Sigma = zeros(size(Y00));
    bilin23_even_SigmaCL = zeros(size(Y00));
    bilin23_even_SigmaD = zeros(size(Y00));
    bilin23_even_M = zeros(size(Y00));
    bilin23_even_MCl = zeros(size(Y00));
    trilin_tau2_even_Sigma = zeros(size(Y00));
    trilin_tau2_even_SigmaCL = zeros(size(Y00));
    trilin_tau2_even_SigmaD = zeros(size(Y00));
    trilin_tau2_even_M = zeros(size(Y00));
    trilin_tau2_even_MCl = zeros(size(Y00));
    trilin_tau3_even_Sigma = zeros(size(Y00));
    trilin_tau3_even_SigmaCL = zeros(size(Y00));
    trilin_tau3_even_SigmaD = zeros(size(Y00));
    trilin_tau3_even_M = zeros(size(Y00));
    trilin_tau3_even_MCl = zeros(size(Y00));
    trilin_tau4_Sigma = zeros(size(Y00));
    trilin_tau4_SigmaCL = zeros(size(Y00));
    trilin_tau4_SigmaD = zeros(size(Y00));
    trilin_tau4_M = zeros(size(Y00));
    trilin_tau4_MCl = zeros(size(Y00));
    bilin12_odd_Sigma = zeros(size(Y00));
    bilin12_odd_SigmaCL = zeros(size(Y00));
    bilin12_odd_SigmaD = zeros(size(Y00));
    bilin12_odd_M = zeros(size(Y00));
    bilin12_odd_MCl = zeros(size(Y00));
    bilin13_odd_Sigma = zeros(size(Y00));
    bilin13_odd_SigmaCL = zeros(size(Y00));
    bilin13_odd_SigmaD = zeros(size(Y00));
    bilin13_odd_M = zeros(size(Y00));
    bilin13_odd_MCl = zeros(size(Y00));
    bilin23_odd_Sigma = zeros(size(Y00));
    bilin23_odd_SigmaCL = zeros(size(Y00));
    bilin23_odd_SigmaD = zeros(size(Y00));
    bilin23_odd_M = zeros(size(Y00));
    bilin23_odd_MCl = zeros(size(Y00));
    trilin_tau2_odd_Sigma = zeros(size(Y00));
    trilin_tau2_odd_SigmaCL = zeros(size(Y00));
    trilin_tau2_odd_SigmaD = zeros(size(Y00));
    trilin_tau2_odd_M = zeros(size(Y00));
    trilin_tau2_odd_MCl = zeros(size(Y00));
    trilin_tau3_odd_Sigma = zeros(size(Y00));
    trilin_tau3_odd_SigmaCL = zeros(size(Y00));
    trilin_tau3_odd_SigmaD = zeros(size(Y00));
    trilin_tau3_odd_M = zeros(size(Y00));
    trilin_tau3_odd_MCl = zeros(size(Y00));
    trilin_tau1_odd_Sigma = zeros(size(Y00));
    trilin_tau1_odd_SigmaCL = zeros(size(Y00));
    trilin_tau1_odd_SigmaD = zeros(size(Y00));
    trilin_tau1_odd_M = zeros(size(Y00));
    trilin_tau1_odd_MCl = zeros(size(Y00));
    bilin12_even_ECl = zeros(size(Y00));
    bilin12_odd_ECl = zeros(size(Y00));
    bilin13_even_ECl = zeros(size(Y00));
    bilin13_odd_ECl = zeros(size(Y00));
    bilin23_even_ECl = zeros(size(Y00));
    bilin23_odd_ECl = zeros(size(Y00));
    trilin_tau1_odd_ECl = zeros(size(Y00));
    trilin_tau2_odd_ECl = zeros(size(Y00));
    trilin_tau2_even_ECl = zeros(size(Y00));
    trilin_tau3_odd_ECl = zeros(size(Y00));
    trilin_tau3_even_ECl = zeros(size(Y00));
    trilin_tau4_even_ECl = zeros(size(Y00));

    % Loop over all pure states of the decomposition of the 3-qubit state
    % to compute the tripartite T, E, and C beads, the Q-beads, and the id bead
    for i = 1:states_123
        norm_state = norm(PureStates123{i});
        weight = norm_state^2;
        rho_pure_123 = (PureStates123{i}/norm_state)*(PureStates123{i}/norm_state)';

        % First calculate all tripartite, Q-beads, and the id bead for Total Correlation function
        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients = {};
        [scaled_coeff] = traverse_and_multiply(LiSAReal{3}, LiSAScalings{3}, scaled_coefficients, rho_pure_123);

        % Select spherical harmonics to the corresponding coefficients
        % Identity bead
        id_bead = scaled_coeff{1}*Y00;
        id_bead_mixed = id_bead_mixed+weight*id_bead;

        % Bead 1, 2, 3
        bead1 = scaled_coeff{2}*Y1n1+scaled_coeff{3}*Y10+scaled_coeff{4}*Y11;
        bead2 = scaled_coeff{5}*Y1n1+scaled_coeff{6}*Y10+scaled_coeff{7}*Y11;
        bead3 = scaled_coeff{8}*Y1n1+scaled_coeff{9}*Y10+scaled_coeff{10}*Y11;

        bead1_mixed = bead1_mixed+weight*bead1;
        bead2_mixed = bead2_mixed+weight*bead2;
        bead3_mixed = bead3_mixed+weight*bead3;

        % Trilinear terms
        % LISA tensors sorted by rank in the structure
        % For tau1: first three odd rank 1 then 5 even terms rank 3

        trilin_tau1_odd = scaled_coeff{38}*Y1n1+scaled_coeff{39}*Y10+...
            scaled_coeff{40}*Y11+scaled_coeff{41}*Y3n3+scaled_coeff{42}*Y3n2+...
            scaled_coeff{43}*Y3n1+scaled_coeff{44}*Y30+ ...
            scaled_coeff{45}*Y31+scaled_coeff{46}*Y32+scaled_coeff{47}*Y33;

        trilin_tau1_odd_mixed = trilin_tau1_odd_mixed+weight*trilin_tau1_odd;

        trilin_tau2_even = scaled_coeff{51}*Y2n2+scaled_coeff{52}*Y2n1+scaled_coeff{53}*Y20+...
            scaled_coeff{54}*Y21+scaled_coeff{55}*Y22;

        trilin_tau2_odd = scaled_coeff{48}*Y1n1+scaled_coeff{49}*Y10+scaled_coeff{50}*Y11;

        trilin_tau2_even_mixed = trilin_tau2_even_mixed+weight*trilin_tau2_even;
        trilin_tau2_odd_mixed = trilin_tau2_odd_mixed+weight*trilin_tau2_odd;

        trilin_tau3_even = scaled_coeff{59}*Y2n2+scaled_coeff{60}*Y2n1+scaled_coeff{61}*Y20+...
            scaled_coeff{62}*Y21+scaled_coeff{63}*Y22;
        trilin_tau3_odd = scaled_coeff{56}*Y1n1+scaled_coeff{57}*Y10+scaled_coeff{58}*Y11;

        trilin_tau3_even_mixed = trilin_tau3_even_mixed+weight*trilin_tau3_even;
        trilin_tau3_odd_mixed = trilin_tau3_odd_mixed+weight*trilin_tau3_odd;

        trilin_tau4 = scaled_coeff{64}*Y00;
        trilin_tau4_mixed = trilin_tau4_mixed+weight*trilin_tau4;

        % Calculate rho connected and rho compound for each pure state
        [rho_connected,rho_compound] = rho_connected_pure(rho_pure_123,dim_123);

        rho_compound_mixed_123 = rho_compound_mixed_123 + rho_compound*weight;
        rho_connected_mixed_123 = rho_connected_mixed_123 + rho_connected*weight;
        % Calculate the tripartite entanglement norm
        [~,~,~,rho_E123] = split_op_tripartite(rho_connected);
        rho_E123 = removeIdLinear(rho_E123);
        % Normalize to 1 (GHZ state has maximum value 1 accordingly)
        EntNorm = EntNorm + weight*getEntanglementNorm(2*real(trace(rho_E123*rho_E123)));

        % Calculate Connected Bilinear Terms E
        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_E = {};
        [scaled_coeff_E] = traverse_and_multiply(LiSAReal{3},LiSAScalings{3},scaled_coefficients_E,rho_connected);

        % Select spherical harmonics to the corresponding coefficients

        % Trilinear terms
        % LISA tensors sorted by rank in the structure
        % For tau1: first three odd rank 1 then 5 even terms rank 3

        trilin_tau1_odd_E = scaled_coeff_E{38}*Y1n1+scaled_coeff_E{39}*Y10+...
            scaled_coeff_E{40}*Y11+scaled_coeff_E{41}*Y3n3+scaled_coeff_E{42}*Y3n2+...
            scaled_coeff_E{43}*Y3n1+scaled_coeff_E{44}*Y30+ ...
            scaled_coeff_E{45}*Y31+scaled_coeff_E{46}*Y32+scaled_coeff_E{47}*Y33;

        trilin_tau1_odd_mixed_E = trilin_tau1_odd_mixed_E+weight*trilin_tau1_odd_E;

        trilin_tau2_odd_E = scaled_coeff_E{48}*Y1n1+scaled_coeff_E{49}*Y10+scaled_coeff_E{50}*Y11;
        trilin_tau2_even_E = scaled_coeff_E{51}*Y2n2+scaled_coeff_E{52}*Y2n1+scaled_coeff_E{53}*Y20+...
            scaled_coeff_E{54}*Y21+scaled_coeff_E{55}*Y22;

        trilin_tau2_odd_mixed_E = trilin_tau2_odd_mixed_E+weight*trilin_tau2_odd_E;
        trilin_tau2_even_mixed_E = trilin_tau2_even_mixed_E+weight*trilin_tau2_even_E;

        trilin_tau3_odd_E = scaled_coeff_E{56}*Y1n1+scaled_coeff_E{57}*Y10+scaled_coeff_E{58}*Y11;
        trilin_tau3_even_E = scaled_coeff_E{59}*Y2n2+scaled_coeff_E{60}*Y2n1+scaled_coeff_E{61}*Y20+...
            scaled_coeff_E{62}*Y21+scaled_coeff_E{63}*Y22;

        trilin_tau3_odd_mixed_E = trilin_tau3_odd_mixed_E+weight*trilin_tau3_odd_E;
        trilin_tau3_even_mixed_E = trilin_tau3_even_mixed_E+weight*trilin_tau3_even_E;

        trilin_tau4_E = scaled_coeff_E{64}*Y00;
        trilin_tau4_mixed_E = trilin_tau4_mixed_E+weight*trilin_tau4_E;

        % Bilinear terms for rho compound C
        % Get scaled coefficients of Lisa tensor basis

        scaled_coefficients_C = {};
        [scaled_coeff_C] = traverse_and_multiply(LiSAReal{3},LiSAScalings{3},scaled_coefficients_C,rho_compound);

        % Trilinear terms
        % LISA tensors sorted by rank in the structure
        % For tau1: first three odd rank 1 then 5 even terms rank 3

        trilin_tau1_odd_C = scaled_coeff_C{38}*Y1n1+scaled_coeff_C{39}*Y10+...
            scaled_coeff_C{40}*Y11+scaled_coeff_C{41}*Y3n3+scaled_coeff_C{42}*Y3n2+...
            scaled_coeff_C{43}*Y3n1+scaled_coeff_C{44}*Y30+ ...
            scaled_coeff_C{45}*Y31+scaled_coeff_C{46}*Y32+scaled_coeff_C{47}*Y33;

        trilin_tau1_odd_mixed_C = trilin_tau1_odd_mixed_C+weight*trilin_tau1_odd_C;

        trilin_tau2_odd_C = scaled_coeff_C{48}*Y1n1+scaled_coeff_C{49}*Y10+scaled_coeff_C{50}*Y11;
        trilin_tau2_even_C = scaled_coeff_C{51}*Y2n2+scaled_coeff_C{52}*Y2n1+scaled_coeff_C{53}*Y20+...
            scaled_coeff_C{54}*Y21+scaled_coeff_C{55}*Y22;

        trilin_tau2_odd_mixed_C = trilin_tau2_odd_mixed_C+weight*trilin_tau2_odd_C;
        trilin_tau2_even_mixed_C = trilin_tau2_even_mixed_C+weight*trilin_tau2_even_C;

        trilin_tau3_odd_C = scaled_coeff_C{56}*Y1n1+scaled_coeff_C{57}*Y10+scaled_coeff_C{58}*Y11;
        trilin_tau3_even_C = scaled_coeff_C{59}*Y2n2+scaled_coeff_C{60}*Y2n1+scaled_coeff_C{61}*Y20+...
            scaled_coeff_C{62}*Y21+scaled_coeff_C{63}*Y22;

        trilin_tau3_odd_mixed_C = trilin_tau3_odd_mixed_C+weight*trilin_tau3_odd_C;
        trilin_tau3_even_mixed_C = trilin_tau3_even_mixed_C+weight*trilin_tau3_even_C;

        trilin_tau4_C = scaled_coeff_C{64}*Y00;
        trilin_tau4_mixed_C = trilin_tau4_mixed_C+weight*trilin_tau4_C;

        rho_original = rho_original+PureStates123{i}*PureStates123{i}';
    end

    % Loops over all pure states of the decompositions of the marginals
    % to compute the bilinear T, E, and C beads of the marginals
    % Marginal of qubits 1 and 2
    for i = 1:states_12
        norm_state = norm(PureStates12{i});
        weight = norm_state^2;
        rho_pure_12 = (PureStates12{i}/norm_state)*(PureStates12{i}/norm_state)';

        % First calculate all bipartite beads for Total Correlation function
        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients = {};
        [scaled_coeff] = traverse_and_multiply(LiSAReal{2}, LiSAScalings{2}, scaled_coefficients, rho_pure_12);

        % Bilinear terms
        % T_00 even and T_2x 5 terms
        bilin12_even = scaled_coeff{8}*Y00 + scaled_coeff{12}*Y2n2+...
            scaled_coeff{13}*Y2n1+scaled_coeff{14}*Y20+...
            scaled_coeff{15}*Y21+scaled_coeff{16}*Y22;

        bilin12_odd = scaled_coeff{9}*Y1n1+scaled_coeff{10}*Y10+scaled_coeff{11}*Y11;

        bilin12_even_mixed = bilin12_even_mixed + weight*bilin12_even;
        bilin12_odd_mixed = bilin12_odd_mixed + weight*bilin12_odd;

        % Calculate rho connected and rho compound for each pure state
        [rho_connected_12, rho_compound_12] = rho_connected_pure(rho_pure_12, dim_12);

        rho_compound_mixed_12 = rho_compound_mixed_12 + rho_compound_12*weight;
        rho_connected_mixed_12 = rho_connected_mixed_12 + rho_connected_12*weight;

        % Bilin 12 Connected
        scaled_coefficients_E = {};
        [scaled_coeff_E] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_E,rho_connected_12);

        bilin12_even_E = scaled_coeff_E{8}*Y00+ scaled_coeff_E{12}*Y2n2+...
            scaled_coeff_E{13}*Y2n1+scaled_coeff_E{14}*Y20+...
            scaled_coeff_E{15}*Y21+scaled_coeff_E{16}*Y22;

        bilin12_odd_E = scaled_coeff_E{9}*Y1n1+scaled_coeff_E{10}*Y10+scaled_coeff_E{11}*Y11;

        bilin12_even_mixed_E = bilin12_even_mixed_E+weight*bilin12_even_E;
        bilin12_odd_mixed_E = bilin12_odd_mixed_E+weight*bilin12_odd_E;

        % Bilin 12 Compound
        scaled_coefficients_C = {};
        [scaled_coeff_C] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_C,rho_compound_12);

        bilin12_even_C = scaled_coeff_C{8}*Y00+ scaled_coeff_C{12}*Y2n2+...
            scaled_coeff_C{13}*Y2n1+scaled_coeff_C{14}*Y20+...
            scaled_coeff_C{15}*Y21+scaled_coeff_C{16}*Y22;

        bilin12_odd_C = scaled_coeff_C{9}*Y1n1+scaled_coeff_C{10}*Y10+scaled_coeff_C{11}*Y11;

        bilin12_even_mixed_C = bilin12_even_mixed_C+weight*bilin12_even_C;
        bilin12_odd_mixed_C = bilin12_odd_mixed_C+weight*bilin12_odd_C;

        rho_12 = rho_12+PureStates12{i}*PureStates12{i}';

    end

    % Marginal of qubits 1 and 3
    for i = 1:states_13
        norm_state = norm(PureStates13{i});
        weight = norm_state^2;
        rho_pure_13 = (PureStates13{i}/norm_state)*(PureStates13{i}/norm_state)';

        % First calculate all bipartite beads for Total Correlation function
        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients = {};
        [scaled_coeff] = traverse_and_multiply(LiSAReal{2}, LiSAScalings{2}, scaled_coefficients, rho_pure_13);

        % Bilinear terms
        % T_00 even and T_2x 5 terms
        bilin13_even = scaled_coeff{8}*Y00+ scaled_coeff{12}*Y2n2+...
            scaled_coeff{13}*Y2n1+scaled_coeff{14}*Y20+...
            scaled_coeff{15}*Y21+scaled_coeff{16}*Y22;

        bilin13_odd = scaled_coeff{9}*Y1n1+scaled_coeff{10}*Y10+scaled_coeff{11}*Y11;

        bilin13_even_mixed = bilin13_even_mixed + weight*bilin13_even;
        bilin13_odd_mixed = bilin13_odd_mixed + weight*bilin13_odd;

        % Calculate rho connected and rho compound for each pure state
        [rho_connected_13, rho_compound_13] = rho_connected_pure(rho_pure_13,dim_13);

        rho_compound_mixed_13 = rho_compound_mixed_13 + rho_compound_13*weight;
        rho_connected_mixed_13 = rho_connected_mixed_13 + rho_connected_13*weight;

        % Bilin 13 Connected
        scaled_coefficients_E = {};
        [scaled_coeff_E] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_E,rho_connected_13);

        bilin13_even_E = scaled_coeff_E{8}*Y00+ scaled_coeff_E{12}*Y2n2+...
            scaled_coeff_E{13}*Y2n1+scaled_coeff_E{14}*Y20+...
            scaled_coeff_E{15}*Y21+scaled_coeff_E{16}*Y22;

        bilin13_odd_E = scaled_coeff_E{9}*Y1n1+scaled_coeff_E{10}*Y10+scaled_coeff_E{11}*Y11;

        bilin13_even_mixed_E = bilin13_even_mixed_E+weight*bilin13_even_E;
        bilin13_odd_mixed_E = bilin13_odd_mixed_E+weight*bilin13_odd_E;

        % Bilin 13 Compound
        scaled_coefficients_C = {};
        [scaled_coeff_C] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_C,rho_compound_13);

        bilin13_even_C = scaled_coeff_C{8}*Y00+ scaled_coeff_C{12}*Y2n2+...
            scaled_coeff_C{13}*Y2n1+scaled_coeff_C{14}*Y20+...
            scaled_coeff_C{15}*Y21+scaled_coeff_C{16}*Y22;

        bilin13_odd_C = scaled_coeff_C{9}*Y1n1+scaled_coeff_C{10}*Y10+scaled_coeff_C{11}*Y11;

        bilin13_even_mixed_C = bilin13_even_mixed_C+weight*bilin13_even_C;
        bilin13_odd_mixed_C = bilin13_odd_mixed_C+weight*bilin13_odd_C;

        rho_13 = rho_13+PureStates13{i}*PureStates13{i}';
    end

    % Marginal of qubits 2 and 3
    for i = 1:states_23
        norm_state = norm(PureStates23{i});
        weight = norm_state^2;
        rho_pure_23 = (PureStates23{i}/norm_state)*(PureStates23{i}/norm_state)';

        % First calculate all bipartite beads for Total Correlation function
        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients = {};
        [scaled_coeff] = traverse_and_multiply(LiSAReal{2}, LiSAScalings{2}, scaled_coefficients, rho_pure_23);

        % Bilin 23
        % T_00 even and T_2x 5 terms
        bilin23_even = scaled_coeff{8}*Y00+ scaled_coeff{12}*Y2n2+...
            scaled_coeff{13}*Y2n1+scaled_coeff{14}*Y20+...
            scaled_coeff{15}*Y21+scaled_coeff{16}*Y22;

        bilin23_odd = scaled_coeff{9}*Y1n1+scaled_coeff{10}*Y10+scaled_coeff{11}*Y11;

        bilin23_even_mixed = bilin23_even_mixed+weight*bilin23_even;
        bilin23_odd_mixed = bilin23_odd_mixed+weight*bilin23_odd;

        % Calculate rho connected and rho compound for each pure state
        [rho_connected_23, rho_compound_23] = rho_connected_pure(rho_pure_23, dim_23);

        rho_compound_mixed_23 = rho_compound_mixed_23 + rho_compound_23*weight;
        rho_connected_mixed_23 = rho_connected_mixed_23 + rho_connected_23*weight;

        % Bilin 23 Connected
        scaled_coefficients_E = {};
        [scaled_coeff_E] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_E,rho_connected_23);

        bilin23_even_E = scaled_coeff_E{8}*Y00+ scaled_coeff_E{12}*Y2n2+...
            scaled_coeff_E{13}*Y2n1+scaled_coeff_E{14}*Y20+...
            scaled_coeff_E{15}*Y21+scaled_coeff_E{16}*Y22;

        bilin23_odd_E = scaled_coeff_E{9}*Y1n1+scaled_coeff_E{10}*Y10+scaled_coeff_E{11}*Y11;

        bilin23_even_mixed_E = bilin23_even_mixed_E+weight*bilin23_even_E;
        bilin23_odd_mixed_E = bilin23_odd_mixed_E+weight*bilin23_odd_E;

        % Bilin 23 Compound
        scaled_coefficients_C = {};
        [scaled_coeff_C] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_C,rho_compound_23);

        bilin23_even_C = scaled_coeff_C{8}*Y00+ scaled_coeff_C{12}*Y2n2+...
            scaled_coeff_C{13}*Y2n1+scaled_coeff_C{14}*Y20+...
            scaled_coeff_C{15}*Y21+scaled_coeff_C{16}*Y22;

        bilin23_odd_C = scaled_coeff_C{9}*Y1n1+scaled_coeff_C{10}*Y10+scaled_coeff_C{11}*Y11;

        bilin23_even_mixed_C = bilin23_even_mixed_C+weight*bilin23_even_C;
        bilin23_odd_mixed_C = bilin23_odd_mixed_C+weight*bilin23_odd_C;

        rho_23 = rho_23+PureStates23{i}*PureStates23{i}';
    end

    % Get two-tangles and three-tangle for connector lines in visualize codes
    two_tangles = zeros(1,3);
    [~,two_tangle12,two_tangle13,two_tangle23,~] = get_three_tangle(rho_original);
    two_tangles(1) = two_tangle12;
    two_tangles(2) = two_tangle13;
    two_tangles(3) = two_tangle23;
    three_tangle = get_three_tangle_mixed(PureStates123);

    % Beads of the separable structure, only required for the corresponding
    % bead types
    if NeedSigma
        % Get rho after subtracting non correlated part (redundant information)
        rho_sigma_123 = compound_analysis(rho_original, rho_compound_mixed_123, rho_connected_mixed_123, dim_123);
        rho_sigma_12 = compound_analysis(rho_12, rho_compound_mixed_12, rho_connected_mixed_12, dim_12);
        rho_sigma_13 = compound_analysis(rho_13, rho_compound_mixed_13, rho_connected_mixed_13, dim_13);
        rho_sigma_23 = compound_analysis(rho_23, rho_compound_mixed_23, rho_connected_mixed_23, dim_23);

        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_Sigma = {};
        [scaled_coeff_Sigma] = traverse_and_multiply(LiSAReal{3},LiSAScalings{3},scaled_coefficients_Sigma,rho_sigma_123);
        scaled_coefficients_Sigma_12 = {};
        [scaled_coeff_Sigma_12] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_Sigma_12,rho_sigma_12);
        scaled_coefficients_Sigma_13 = {};
        [scaled_coeff_Sigma_13] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_Sigma_13,rho_sigma_13);
        scaled_coefficients_Sigma_23 = {};
        [scaled_coeff_Sigma_23] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_Sigma_23,rho_sigma_23);

        % Select spherical harmonics to the corresponding coefficients
        % Bilin 12
        % T_00 even and T_2x 5 terms
        bilin12_even_Sigma = scaled_coeff_Sigma_12{8}*Y00 + scaled_coeff_Sigma_12{12}*Y2n2+...
            scaled_coeff_Sigma_12{13}*Y2n1 + scaled_coeff_Sigma_12{14}*Y20+...
            scaled_coeff_Sigma_12{15}*Y21 + scaled_coeff_Sigma_12{16}*Y22;

        bilin12_odd_Sigma = scaled_coeff_Sigma_12{9}*Y1n1+scaled_coeff_Sigma_12{10}*Y10+scaled_coeff_Sigma_12{11}*Y11;

        % Bilin 13
        % T_00 even and T_2x 5 terms
        bilin13_even_Sigma = scaled_coeff_Sigma_13{8}*Y00+ scaled_coeff_Sigma_13{12}*Y2n2+...
            scaled_coeff_Sigma_13{13}*Y2n1+scaled_coeff_Sigma_13{14}*Y20+...
            scaled_coeff_Sigma_13{15}*Y21+scaled_coeff_Sigma_13{16}*Y22;

        bilin13_odd_Sigma = scaled_coeff_Sigma_13{9}*Y1n1+scaled_coeff_Sigma_13{10}*Y10+scaled_coeff_Sigma_13{11}*Y11;

        % Bilin 23
        % T_00 even and T_2x 5 terms
        bilin23_even_Sigma = scaled_coeff_Sigma_23{8}*Y00 + scaled_coeff_Sigma_23{12}*Y2n2 +...
            scaled_coeff_Sigma_23{13}*Y2n1 + scaled_coeff_Sigma_23{14}*Y20 +...
            scaled_coeff_Sigma_23{15}*Y21 + scaled_coeff_Sigma_23{16}*Y22;

        bilin23_odd_Sigma = scaled_coeff_Sigma_23{9}*Y1n1+scaled_coeff_Sigma_23{10}*Y10+scaled_coeff_Sigma_23{11}*Y11;

        % Trilinear terms
        % LISA tensors sorted by rank in the structure
        % For tau1: first three odd rank 1 then 5 even terms rank 3

        trilin_tau1_odd_Sigma = scaled_coeff_Sigma{38}*Y1n1+scaled_coeff_Sigma{39}*Y10+...
            scaled_coeff_Sigma{40}*Y11+scaled_coeff_Sigma{41}*Y3n3+scaled_coeff_Sigma{42}*Y3n2+...
            scaled_coeff_Sigma{43}*Y3n1+scaled_coeff_Sigma{44}*Y30+ ...
            scaled_coeff_Sigma{45}*Y31+scaled_coeff_Sigma{46}*Y32+scaled_coeff_Sigma{47}*Y33;

        trilin_tau2_odd_Sigma = scaled_coeff_Sigma{48}*Y1n1+scaled_coeff_Sigma{49}*Y10+scaled_coeff_Sigma{50}*Y11;
        trilin_tau2_even_Sigma = scaled_coeff_Sigma{51}*Y2n2+scaled_coeff_Sigma{52}*Y2n1+scaled_coeff_Sigma{53}*Y20+...
            scaled_coeff_Sigma{54}*Y21+scaled_coeff_Sigma{55}*Y22;

        trilin_tau3_odd_Sigma = scaled_coeff_Sigma{56}*Y1n1+scaled_coeff_Sigma{57}*Y10+scaled_coeff_Sigma{58}*Y11;
        trilin_tau3_even_Sigma = scaled_coeff_Sigma{59}*Y2n2+scaled_coeff_Sigma{60}*Y2n1+scaled_coeff_Sigma{61}*Y20+...
            scaled_coeff_Sigma{62}*Y21+scaled_coeff_Sigma{63}*Y22;

        trilin_tau4_Sigma = scaled_coeff_Sigma{64}*Y00;
    end

    % Mutual structure and its beads, only required for the corresponding
    % bead types
    if NeedM
        rho_M = rho_mutualinformation(rho_original);
        rho_M12 = rho_mutualinformation(rho_12);
        rho_M13 = rho_mutualinformation(rho_13);
        rho_M23 = rho_mutualinformation(rho_23);

        % M-beads: total mutual information
        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_M = {};
        [scaled_coeff_M] = traverse_and_multiply(LiSAReal{3},LiSAScalings{3},scaled_coefficients_M,rho_M);

        % Select spherical harmonics to the corresponding coefficients

        % Bilin 12
        % T_00 even and T_2x 5 terms
        bilin12_even_M= scaled_coeff_M{11}*Y00+ scaled_coeff_M{15}*Y2n2+...
            scaled_coeff_M{16}*Y2n1+scaled_coeff_M{17}*Y20+...
            scaled_coeff_M{18}*Y21+scaled_coeff_M{19}*Y22;

        bilin12_odd_M = scaled_coeff_M{12}*Y1n1+scaled_coeff_M{13}*Y10+scaled_coeff_M{14}*Y11;

        % Bilin 13
        % T_00 even and T_2x 5 terms
        bilin13_even_M= scaled_coeff_M{20}*Y00+ scaled_coeff_M{24}*Y2n2+...
            scaled_coeff_M{25}*Y2n1+scaled_coeff_M{26}*Y20+...
            scaled_coeff_M{27}*Y21+scaled_coeff_M{28}*Y22;

        bilin13_odd_M = scaled_coeff_M{21}*Y1n1+scaled_coeff_M{22}*Y10+scaled_coeff_M{23}*Y11;

        % Bilin 23
        % T_00 even and T_2x 5 terms
        bilin23_even_M= scaled_coeff_M{29}*Y00+ scaled_coeff_M{33}*Y2n2+...
            scaled_coeff_M{34}*Y2n1+scaled_coeff_M{35}*Y20+...
            scaled_coeff_M{36}*Y21+scaled_coeff_M{37}*Y22;

        bilin23_odd_M = scaled_coeff_M{30}*Y1n1+scaled_coeff_M{31}*Y10+scaled_coeff_M{32}*Y11;

        % Trilinear terms
        % LISA tensors sorted by rank in the structure
        % For tau1: first three odd rank 1 then 5 even terms rank 3

        trilin_tau1_odd_M = scaled_coeff_M{38}*Y1n1+scaled_coeff_M{39}*Y10+...
            scaled_coeff_M{40}*Y11+scaled_coeff_M{41}*Y3n3+scaled_coeff_M{42}*Y3n2+...
            scaled_coeff_M{43}*Y3n1+scaled_coeff_M{44}*Y30+ ...
            scaled_coeff_M{45}*Y31+scaled_coeff_M{46}*Y32+scaled_coeff_M{47}*Y33;

        trilin_tau2_odd_M = scaled_coeff_M{48}*Y1n1+scaled_coeff_M{49}*Y10+scaled_coeff_M{50}*Y11;
        trilin_tau2_even_M = scaled_coeff_M{51}*Y2n2+scaled_coeff_M{52}*Y2n1+scaled_coeff_M{53}*Y20+...
            scaled_coeff_M{54}*Y21+scaled_coeff_M{55}*Y22;

        trilin_tau3_odd_M = scaled_coeff_M{56}*Y1n1+scaled_coeff_M{57}*Y10+scaled_coeff_M{58}*Y11;
        trilin_tau3_even_M = scaled_coeff_M{59}*Y2n2+scaled_coeff_M{60}*Y2n1+scaled_coeff_M{61}*Y20+...
            scaled_coeff_M{62}*Y21+scaled_coeff_M{63}*Y22;

        trilin_tau4_M = scaled_coeff_M{64}*Y00;
    end

    % Measurement projectors, classical beads, nonclassical beads, only
    % required for the corresponding bead types
    if NeedProj
        % Split up the separable and the mutual structure into 3 bilinear and
        % 1 trilinear part to independently find the projectors for the
        % classical correlation
        [~,~,~,rho_sigma123] = split_op_tripartite(rho_sigma_123);
        [~,~,~,rho_M123] = split_op_tripartite(rho_M);

        % Check each part separately, the projectors are determined
        % recursively for the marginals and the tripartite structure
        % (OptStructure only affects the tripartite structure)
        if ClassCorrRecursive
            % Marginal of qubits 1 and 2, the target qubit within the marginal
            % is the target qubit itself
            if TargetQubits(1) == 1
                targetqubit_12 = 1;
            elseif TargetQubits(1) == 2
                targetqubit_12 = 2;
            end

            [rho_sigmaCl12,rho_MCl12, ~, rho_ECl12] = find_projectors_bipartite(rho_12, rho_sigma_12, rho_M12, rho_connected_mixed_12, targetqubit_12, OptStructure, MeasurementDirection);

            % Marginal of qubits 1 and 3, the target qubit can only be 1 or 3. If
            % it is 3 it has to be reindexed with 2 since the marginal is a
            % two-qubit state
            if TargetQubits(2) == 3
                targetqubit_13 = 2;
            elseif TargetQubits(2) == 1
                targetqubit_13 = 1;
            end

            [rho_sigmaCl13,rho_MCl13, ~, rho_ECl13] = find_projectors_bipartite(rho_13, rho_sigma_13, rho_M13, rho_connected_mixed_13, targetqubit_13, OptStructure, MeasurementDirection);

            % Marginal of qubits 2 and 3, the target qubit can only be 2 or 3. If
            % it is 3 it has to be reindexed with 2 and if it is 2 with 1 since
            % the marginal is a two-qubit state
            if TargetQubits(3) == 3
                targetqubit_23 = 2;
            elseif TargetQubits(3) == 2
                targetqubit_23 = 1;
            end

            [rho_sigmaCl23,rho_MCl23, ~, rho_ECl23] = find_projectors_bipartite(rho_23,rho_sigma_23,rho_M23,rho_connected_mixed_23, targetqubit_23, OptStructure, MeasurementDirection);

            % Isolate genuine tripartite correlations, remove any other components
            % of the mutual structure
            [~,~,~,rho_123_tri] = split_op_tripartite(rho_M);
            rho_123_tri = removeIdLinear(rho_123_tri, false);

            % Find the projectors for the tripartite structure for both target
            % qubits
            [rho_sigmaCl123,rho_MCl123,rho_ECl123] = find_projectors_tripartite(rho_123_tri,rho_sigma123,rho_M123,rho_connected_mixed_123,TargetQubits(4:5), OptStructure);
        else
            % Determine the projectors for the full structure and apply them to
            % the marginals as well
            [rho_post, rho_post_ent, rho_post_sep] = opt_classcor_universal3(rho_original, TargetQubits(4:5), ...
                                                                            rho_M123, rho_M12, ...
                                                                            rho_M13, rho_M23, ...
                                                                            rho_connected_mixed_123, rho_connected_mixed_12, ...
                                                                            rho_connected_mixed_13, rho_connected_mixed_23, ...
                                                                            rho_sigma123, rho_sigma_12, ...
                                                                            rho_sigma_13, rho_sigma_23);

            % Extract the classical structures from the post-measurement structures
            rho_MCl123 = rho_post{1};
            rho_MCl12 = rho_post{2};
            rho_MCl13 = rho_post{3};
            rho_MCl23 = rho_post{4};

            rho_ECl123 = rho_post_ent{1};
            rho_ECl12 = rho_post_ent{2};
            rho_ECl13 = rho_post_ent{3};
            rho_ECl23 = rho_post_ent{4};

            rho_sigmaCl123 = rho_post_sep{1};
            rho_sigmaCl12 = rho_post_sep{2};
            rho_sigmaCl13 = rho_post_sep{3};
            rho_sigmaCl23 = rho_post_sep{4};
        end

        % Nonclassical part of the separable structure
        rho_sigmaD_123 = rho_sigma_123 - rho_sigmaCl123;
        rho_sigmaD_12 = rho_sigma_12 - rho_sigmaCl12;
        rho_sigmaD_13 = rho_sigma_13 - rho_sigmaCl13;
        rho_sigmaD_23 = rho_sigma_23 - rho_sigmaCl23;

        % Beads of the classical part of the connected structure
        scaled_coefficients_ECl = {};
        [scaled_coeff_ECl] = traverse_and_multiply(LiSAReal{3},LiSAScalings{3},scaled_coefficients_ECl,rho_ECl123);
        scaled_coefficients_ECl_12 = {};
        [scaled_coeff_ECl_12] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_ECl_12,rho_ECl12);
        scaled_coefficients_ECl_13 = {};
        [scaled_coeff_ECl_13] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_ECl_13,rho_ECl13);
        scaled_coefficients_ECl_23 = {};
        [scaled_coeff_ECl_23] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_ECl_23,rho_ECl23);

        bilin12_even_ECl = scaled_coeff_ECl_12{8}*Y00+ scaled_coeff_ECl_12{12}*Y2n2+...
            scaled_coeff_ECl_12{13}*Y2n1+scaled_coeff_ECl_12{14}*Y20+...
            scaled_coeff_ECl_12{15}*Y21+scaled_coeff_ECl_12{16}*Y22;

        bilin12_odd_ECl = scaled_coeff_ECl_12{9}*Y1n1+scaled_coeff_ECl_12{10}*Y10+scaled_coeff_ECl_12{11}*Y11;

        bilin13_even_ECl = scaled_coeff_ECl_13{8}*Y00+ scaled_coeff_ECl_13{12}*Y2n2+...
            scaled_coeff_ECl_13{13}*Y2n1+scaled_coeff_ECl_13{14}*Y20+...
            scaled_coeff_ECl_13{15}*Y21+scaled_coeff_ECl_13{16}*Y22;

        bilin13_odd_ECl = scaled_coeff_ECl_13{9}*Y1n1+scaled_coeff_ECl_13{10}*Y10+scaled_coeff_ECl_13{11}*Y11;

        bilin23_even_ECl = scaled_coeff_ECl_23{8}*Y00+ scaled_coeff_ECl_23{12}*Y2n2+...
            scaled_coeff_ECl_23{13}*Y2n1+scaled_coeff_ECl_23{14}*Y20+...
            scaled_coeff_ECl_23{15}*Y21+scaled_coeff_ECl_23{16}*Y22;

        bilin23_odd_ECl = scaled_coeff_ECl_23{9}*Y1n1+scaled_coeff_ECl_23{10}*Y10+scaled_coeff_ECl_23{11}*Y11;

        trilin_tau1_odd_ECl = scaled_coeff_ECl{38}*Y1n1+scaled_coeff_ECl{39}*Y10+...
            scaled_coeff_ECl{40}*Y11+scaled_coeff_ECl{41}*Y3n3+scaled_coeff_ECl{42}*Y3n2+...
            scaled_coeff_ECl{43}*Y3n1+scaled_coeff_ECl{44}*Y30+ ...
            scaled_coeff_ECl{45}*Y31+scaled_coeff_ECl{46}*Y32+scaled_coeff_ECl{47}*Y33;

        trilin_tau2_odd_ECl = scaled_coeff_ECl{48}*Y1n1+scaled_coeff_ECl{49}*Y10+scaled_coeff_ECl{50}*Y11;
        trilin_tau2_even_ECl = scaled_coeff_ECl{51}*Y2n2+scaled_coeff_ECl{52}*Y2n1+scaled_coeff_ECl{53}*Y20+...
            scaled_coeff_ECl{54}*Y21+scaled_coeff_ECl{55}*Y22;

        trilin_tau3_odd_ECl = scaled_coeff_ECl{56}*Y1n1+scaled_coeff_ECl{57}*Y10+scaled_coeff_ECl{58}*Y11;
        trilin_tau3_even_ECl = scaled_coeff_ECl{59}*Y2n2+scaled_coeff_ECl{60}*Y2n1+scaled_coeff_ECl{61}*Y20+...
            scaled_coeff_ECl{62}*Y21+scaled_coeff_ECl{63}*Y22;

        trilin_tau4_even_ECl = scaled_coeff_ECl{64}*Y00;

        % Beads of the classical part of the separable structure
        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_SigmaCL = {};
        [scaled_coeff_SigmaCL] = traverse_and_multiply(LiSAReal{3},LiSAScalings{3},scaled_coefficients_SigmaCL,rho_sigmaCl123);
        scaled_coefficients_SigmaCL_12 = {};
        [scaled_coeff_SigmaCL_12] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_SigmaCL_12,rho_sigmaCl12);
        scaled_coefficients_SigmaCL_13 = {};
        [scaled_coeff_SigmaCL_13] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_SigmaCL_13,rho_sigmaCl13);
        scaled_coefficients_SigmaCL_23 = {};
        [scaled_coeff_SigmaCL_23] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_SigmaCL_23,rho_sigmaCl23);

        % Select spherical harmonics to the corresponding coefficients

        % Bilin 12
        % T_00 even and T_2x 5 terms
        bilin12_even_SigmaCL = scaled_coeff_SigmaCL_12{8}*Y00+ scaled_coeff_SigmaCL_12{12}*Y2n2+...
            scaled_coeff_SigmaCL_12{13}*Y2n1+scaled_coeff_SigmaCL_12{14}*Y20+...
            scaled_coeff_SigmaCL_12{15}*Y21+scaled_coeff_SigmaCL_12{16}*Y22;

        bilin12_odd_SigmaCL = scaled_coeff_SigmaCL_12{9}*Y1n1+scaled_coeff_SigmaCL_12{10}*Y10+scaled_coeff_SigmaCL_12{11}*Y11;

        % Bilin 13
        % T_00 even and T_2x 5 terms
        bilin13_even_SigmaCL = scaled_coeff_SigmaCL_13{8}*Y00+ scaled_coeff_SigmaCL_13{12}*Y2n2+...
            scaled_coeff_SigmaCL_13{13}*Y2n1+scaled_coeff_SigmaCL_13{14}*Y20+...
            scaled_coeff_SigmaCL_13{15}*Y21+scaled_coeff_SigmaCL_13{16}*Y22;

        bilin13_odd_SigmaCL = scaled_coeff_SigmaCL_13{9}*Y1n1+scaled_coeff_SigmaCL_13{10}*Y10+scaled_coeff_SigmaCL_13{11}*Y11;

        % Bilin 23
        % T_00 even and T_2x 5 terms
        bilin23_even_SigmaCL = scaled_coeff_SigmaCL_23{8}*Y00+ scaled_coeff_SigmaCL_23{12}*Y2n2+...
            scaled_coeff_SigmaCL_23{13}*Y2n1+scaled_coeff_SigmaCL_23{14}*Y20+...
            scaled_coeff_SigmaCL_23{15}*Y21+scaled_coeff_SigmaCL_23{16}*Y22;

        bilin23_odd_SigmaCL = scaled_coeff_SigmaCL_23{9}*Y1n1+scaled_coeff_SigmaCL_23{10}*Y10+scaled_coeff_SigmaCL_23{11}*Y11;

        % Trilinear terms
        % LISA tensors sorted by rank in the structure
        % For tau1: first three odd rank 1 then 5 even terms rank 3

        trilin_tau1_odd_SigmaCL = scaled_coeff_SigmaCL{38}*Y1n1+scaled_coeff_SigmaCL{39}*Y10+...
            scaled_coeff_SigmaCL{40}*Y11+scaled_coeff_SigmaCL{41}*Y3n3+scaled_coeff_SigmaCL{42}*Y3n2+...
            scaled_coeff_SigmaCL{43}*Y3n1+scaled_coeff_SigmaCL{44}*Y30+ ...
            scaled_coeff_SigmaCL{45}*Y31+scaled_coeff_SigmaCL{46}*Y32+scaled_coeff_SigmaCL{47}*Y33;

        trilin_tau2_odd_SigmaCL = scaled_coeff_SigmaCL{48}*Y1n1+scaled_coeff_SigmaCL{49}*Y10+scaled_coeff_SigmaCL{50}*Y11;
        trilin_tau2_even_SigmaCL = scaled_coeff_SigmaCL{51}*Y2n2+scaled_coeff_SigmaCL{52}*Y2n1+scaled_coeff_SigmaCL{53}*Y20+...
            scaled_coeff_SigmaCL{54}*Y21+scaled_coeff_SigmaCL{55}*Y22;

        trilin_tau3_odd_SigmaCL = scaled_coeff_SigmaCL{56}*Y1n1+scaled_coeff_SigmaCL{57}*Y10+scaled_coeff_SigmaCL{58}*Y11;
        trilin_tau3_even_SigmaCL = scaled_coeff_SigmaCL{59}*Y2n2+scaled_coeff_SigmaCL{60}*Y2n1+scaled_coeff_SigmaCL{61}*Y20+...
            scaled_coeff_SigmaCL{62}*Y21+scaled_coeff_SigmaCL{63}*Y22;

        trilin_tau4_SigmaCL = scaled_coeff_SigmaCL{64}*Y00;

        % Beads of the nonclassical part of the separable structure
        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_SigmaD = {};
        [scaled_coeff_SigmaD] = traverse_and_multiply(LiSAReal{3},LiSAScalings{3},scaled_coefficients_SigmaD, rho_sigmaD_123);
        scaled_coefficients_SigmaD_12 = {};
        [scaled_coeff_SigmaD_12] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_SigmaD_12, rho_sigmaD_12);
        scaled_coefficients_SigmaD_13 = {};
        [scaled_coeff_SigmaD_13] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_SigmaD_13, rho_sigmaD_13);
        scaled_coefficients_SigmaD_23 = {};
        [scaled_coeff_SigmaD_23] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_SigmaD_23, rho_sigmaD_23);

        % Select spherical harmonics to the corresponding coefficients
        % Bilin 12
        % T_00 even and T_2x 5 terms
        bilin12_even_SigmaD = scaled_coeff_SigmaD_12{8}*Y00+ scaled_coeff_SigmaD_12{12}*Y2n2+...
            scaled_coeff_SigmaD_12{13}*Y2n1+scaled_coeff_SigmaD_12{14}*Y20+...
            scaled_coeff_SigmaD_12{15}*Y21+scaled_coeff_SigmaD_12{16}*Y22;

        bilin12_odd_SigmaD = scaled_coeff_SigmaD_12{9}*Y1n1+scaled_coeff_SigmaD_12{10}*Y10+scaled_coeff_SigmaD_12{11}*Y11;

        % Bilin 13
        % T_00 even and T_2x 5 terms
        bilin13_even_SigmaD = scaled_coeff_SigmaD_13{8}*Y00+ scaled_coeff_SigmaD_13{12}*Y2n2+...
            scaled_coeff_SigmaD_13{13}*Y2n1+scaled_coeff_SigmaD_13{14}*Y20+...
            scaled_coeff_SigmaD_13{15}*Y21+scaled_coeff_SigmaD_13{16}*Y22;

        bilin13_odd_SigmaD = scaled_coeff_SigmaD_13{9}*Y1n1+scaled_coeff_SigmaD_13{10}*Y10+scaled_coeff_SigmaD_13{11}*Y11;

        % Bilin 23
        % T_00 even and T_2x 5 terms
        bilin23_even_SigmaD = scaled_coeff_SigmaD_23{8}*Y00+ scaled_coeff_SigmaD_23{12}*Y2n2+...
            scaled_coeff_SigmaD_23{13}*Y2n1+scaled_coeff_SigmaD_23{14}*Y20+...
            scaled_coeff_SigmaD_23{15}*Y21+scaled_coeff_SigmaD_23{16}*Y22;

        bilin23_odd_SigmaD = scaled_coeff_SigmaD_23{9}*Y1n1+scaled_coeff_SigmaD_23{10}*Y10+scaled_coeff_SigmaD_23{11}*Y11;

        % Trilinear terms
        % LISA tensors sorted by rank in the structure
        % For tau1: first three odd rank 1 then 5 even terms rank 3

        trilin_tau1_odd_SigmaD = scaled_coeff_SigmaD{38}*Y1n1+scaled_coeff_SigmaD{39}*Y10+...
            scaled_coeff_SigmaD{40}*Y11+scaled_coeff_SigmaD{41}*Y3n3+scaled_coeff_SigmaD{42}*Y3n2+...
            scaled_coeff_SigmaD{43}*Y3n1+scaled_coeff_SigmaD{44}*Y30+ ...
            scaled_coeff_SigmaD{45}*Y31+scaled_coeff_SigmaD{46}*Y32+scaled_coeff_SigmaD{47}*Y33;

        trilin_tau2_odd_SigmaD = scaled_coeff_SigmaD{48}*Y1n1+scaled_coeff_SigmaD{49}*Y10+scaled_coeff_SigmaD{50}*Y11;
        trilin_tau2_even_SigmaD = scaled_coeff_SigmaD{51}*Y2n2+scaled_coeff_SigmaD{52}*Y2n1+scaled_coeff_SigmaD{53}*Y20+...
            scaled_coeff_SigmaD{54}*Y21+scaled_coeff_SigmaD{55}*Y22;

        trilin_tau3_odd_SigmaD = scaled_coeff_SigmaD{56}*Y1n1+scaled_coeff_SigmaD{57}*Y10+scaled_coeff_SigmaD{58}*Y11;
        trilin_tau3_even_SigmaD = scaled_coeff_SigmaD{59}*Y2n2+scaled_coeff_SigmaD{60}*Y2n1+scaled_coeff_SigmaD{61}*Y20+...
            scaled_coeff_SigmaD{62}*Y21+scaled_coeff_SigmaD{63}*Y22;

        trilin_tau4_SigmaD = scaled_coeff_SigmaD{64}*Y00;

        % Beads of the classical mutual structure
        % Get coefficients for classical mutual information
        scaled_coefficients_MCl = {};
        [scaled_coeff_MCl] = traverse_and_multiply(LiSAReal{3},LiSAScalings{3},scaled_coefficients_MCl, rho_MCl123);
        scaled_coefficients_MCl_12 = {};
        [scaled_coeff_MCl_12] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_MCl_12, rho_MCl12);
        scaled_coefficients_MCl_13 = {};
        [scaled_coeff_MCl_13] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_MCl_13, rho_MCl13);
        scaled_coefficients_MCl_23 = {};
        [scaled_coeff_MCl_23] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_MCl_23, rho_MCl23);

        % Select spherical harmonics to the corresponding coefficients
        % Bilin 12
        % T_00 even and T_2x 5 terms
        bilin12_even_MCl = scaled_coeff_MCl_12{8}*Y00+ scaled_coeff_MCl_12{12}*Y2n2+...
            scaled_coeff_MCl_12{13}*Y2n1+scaled_coeff_MCl_12{14}*Y20+...
            scaled_coeff_MCl_12{15}*Y21+scaled_coeff_MCl_12{16}*Y22;

        bilin12_odd_MCl = scaled_coeff_MCl_12{9}*Y1n1+scaled_coeff_MCl_12{10}*Y10+scaled_coeff_MCl_12{11}*Y11;

        % Bilin 13
        % T_00 even and T_2x 5 terms
        bilin13_even_MCl = scaled_coeff_MCl_13{8}*Y00+ scaled_coeff_MCl_13{12}*Y2n2+...
            scaled_coeff_MCl_13{13}*Y2n1+scaled_coeff_MCl_13{14}*Y20+...
            scaled_coeff_MCl_13{15}*Y21+scaled_coeff_MCl_13{16}*Y22;

        bilin13_odd_MCl = scaled_coeff_MCl_13{9}*Y1n1+scaled_coeff_MCl_13{10}*Y10+scaled_coeff_MCl_13{11}*Y11;

        % Bilin 23
        % T_00 even and T_2x 5 terms
        bilin23_even_MCl = scaled_coeff_MCl_23{8}*Y00+ scaled_coeff_MCl_23{12}*Y2n2+...
            scaled_coeff_MCl_23{13}*Y2n1+scaled_coeff_MCl_23{14}*Y20+...
            scaled_coeff_MCl_23{15}*Y21+scaled_coeff_MCl_23{16}*Y22;

        bilin23_odd_MCl = scaled_coeff_MCl_23{9}*Y1n1+scaled_coeff_MCl_23{10}*Y10+scaled_coeff_MCl_23{11}*Y11;

        % Trilinear terms
        % LISA tensors sorted by rank in the structure
        % For tau1: first three odd rank 1 then 5 even terms rank 3

        trilin_tau1_odd_MCl = scaled_coeff_MCl{38}*Y1n1+scaled_coeff_MCl{39}*Y10+...
            scaled_coeff_MCl{40}*Y11+scaled_coeff_MCl{41}*Y3n3+scaled_coeff_MCl{42}*Y3n2+...
            scaled_coeff_MCl{43}*Y3n1+scaled_coeff_MCl{44}*Y30+ ...
            scaled_coeff_MCl{45}*Y31+scaled_coeff_MCl{46}*Y32+scaled_coeff_MCl{47}*Y33;

        trilin_tau2_odd_MCl = scaled_coeff_MCl{48}*Y1n1+scaled_coeff_MCl{49}*Y10+scaled_coeff_MCl{50}*Y11;
        trilin_tau2_even_MCl = scaled_coeff_MCl{51}*Y2n2+scaled_coeff_MCl{52}*Y2n1+scaled_coeff_MCl{53}*Y20+...
            scaled_coeff_MCl{54}*Y21+scaled_coeff_MCl{55}*Y22;

        trilin_tau3_odd_MCl = scaled_coeff_MCl{56}*Y1n1+scaled_coeff_MCl{57}*Y10+scaled_coeff_MCl{58}*Y11;
        trilin_tau3_even_MCl = scaled_coeff_MCl{59}*Y2n2+scaled_coeff_MCl{60}*Y2n1+scaled_coeff_MCl{61}*Y20+...
            scaled_coeff_MCl{62}*Y21+scaled_coeff_MCl{63}*Y22;

        trilin_tau4_MCl = scaled_coeff_MCl{64}*Y00;
    end

    % Collect beads for output according to EVEN and then ODD case
    Beads = {id_bead_mixed, bead1_mixed, bead2_mixed, bead3_mixed, ...
        bilin12_even_mixed,     bilin12_even_mixed_E,     bilin12_even_mixed_C, bilin12_even_Sigma, ...
        bilin12_even_SigmaCL,        bilin12_even_SigmaD,           bilin12_even_M,       bilin12_even_MCl, ...
        bilin13_even_mixed,     bilin13_even_mixed_E,     bilin13_even_mixed_C, bilin13_even_Sigma, ...
        bilin13_even_SigmaCL,        bilin13_even_SigmaD,           bilin13_even_M,       bilin13_even_MCl, ...
        bilin23_even_mixed,     bilin23_even_mixed_E,     bilin23_even_mixed_C, bilin23_even_Sigma, ...
        bilin23_even_SigmaCL,        bilin23_even_SigmaD,           bilin23_even_M,       bilin23_even_MCl, ...
        trilin_tau2_even_mixed, trilin_tau2_even_mixed_E, trilin_tau2_even_mixed_C, ...
        trilin_tau2_even_Sigma,    trilin_tau2_even_SigmaCL,      trilin_tau2_even_SigmaD, ...
        trilin_tau2_even_M,     trilin_tau2_even_MCl,...
        trilin_tau3_even_mixed, trilin_tau3_even_mixed_E, trilin_tau3_even_mixed_C, ...
        trilin_tau3_even_Sigma,    trilin_tau3_even_SigmaCL,      trilin_tau3_even_SigmaD, ...
        trilin_tau3_even_M,     trilin_tau3_even_MCl, ...
        trilin_tau4_mixed,      trilin_tau4_mixed_E,      trilin_tau4_mixed_C,...
        trilin_tau4_Sigma,trilin_tau4_SigmaCL,trilin_tau4_SigmaD,trilin_tau4_M,trilin_tau4_MCl, ...
        bilin12_odd_mixed, bilin12_odd_mixed_E, bilin12_odd_mixed_C, bilin12_odd_Sigma,bilin12_odd_SigmaCL,bilin12_odd_SigmaD,bilin12_odd_M,bilin12_odd_MCl,...
        bilin13_odd_mixed, bilin13_odd_mixed_E, bilin13_odd_mixed_C,bilin13_odd_Sigma,bilin13_odd_SigmaCL,bilin13_odd_SigmaD,bilin13_odd_M,bilin13_odd_MCl,...
        bilin23_odd_mixed, bilin23_odd_mixed_E,bilin23_odd_mixed_C,bilin23_odd_Sigma,bilin23_odd_SigmaCL,bilin23_odd_SigmaD,bilin23_odd_M,bilin23_odd_MCl,...
        trilin_tau2_odd_mixed, trilin_tau2_odd_mixed_E, trilin_tau2_odd_mixed_C,trilin_tau2_odd_Sigma,trilin_tau2_odd_SigmaCL,trilin_tau2_odd_SigmaD,trilin_tau2_odd_M,trilin_tau2_odd_MCl,...
        trilin_tau3_odd_mixed, trilin_tau3_odd_mixed_E, trilin_tau3_odd_mixed_C,trilin_tau3_odd_Sigma,trilin_tau3_odd_SigmaCL, trilin_tau3_odd_SigmaD,trilin_tau3_odd_M,trilin_tau3_odd_MCl, ...
        trilin_tau1_odd_mixed,trilin_tau1_odd_mixed_E,trilin_tau1_odd_mixed_C,trilin_tau1_odd_Sigma,trilin_tau1_odd_SigmaCL,trilin_tau1_odd_SigmaD,trilin_tau1_odd_M,trilin_tau1_odd_MCl,...
        bilin12_even_ECl, bilin12_odd_ECl, bilin13_even_ECl, bilin13_odd_ECl, bilin23_even_ECl, bilin23_odd_ECl, trilin_tau1_odd_ECl, trilin_tau2_odd_ECl, trilin_tau2_even_ECl, ...
        trilin_tau3_odd_ECl, trilin_tau3_even_ECl, trilin_tau4_even_ECl, TargetQubits,two_tangles,three_tangle};
end
