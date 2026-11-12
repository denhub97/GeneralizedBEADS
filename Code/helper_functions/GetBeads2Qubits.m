% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [Beads, EntNorm] = GetBeads2Qubits(PureStates, BeadTypes, TargetQubit, GridPoints, OptStructure, MeasurementDirection)

    %!!! Compute the spherical function values of all beads
    %    of a 2-qubit state, given as a pure state decomposition
    %    of the density operator. The beads of the pure states are weighted
    %    with their squared norm and summed. Next to the beads of the total
    %    correlation (T), the connected correlations (E) and the compound
    %    part (C), the beads of the separable structure + its classical
    %    and nonclassical part, of the mutual structure and its classical
    %    part, and of the classical part of the connected structure are
    %    computed (using the measurement projectors found for the target
    %    qubit), together with the difference beads of the local states. The
    %    entanglement norm of the state is computed as well
    %
    %    input: PureStates: Cell array of the pure states constituting the
    %                       decomposition (the weights are contained in
    %                       the state vectors, i.e., sqrt(p)*|state>)
    %           BeadTypes: List of the bead types (strings or character
    %                      arrays, at least one element)
    %           TargetQubit: The qubit which is measured for the classical
    %                        structures (1 or 2)
    %           GridPoints: The beads resolution [theta, phi]
    %           OptStructure: Structure which the classical beads are
    %                       optimized for (integer): -1 = entangling
    %                       structure (the overlap with the separable
    %                       structure is minimized), 0 = no structure
    %                       optimization, 1 = separable structure (the
    %                       overlap with the separable structure is
    %                       maximized) (default 0)
    %           MeasurementDirection: Measurement direction [x y z] used
    %                       by find_projectors_bipartite if the optimal
    %                       orientation is arbitrary. If empty (default),
    %                       the direction is requested in the command window
    %
    %    output: Beads: Cell array with the beads, in the order: identity
    %                   bead, Q-beads 1 and 2, the bilinear even
    %                   beads (T, E, C, separable structure, the classical
    %                   and nonclassical separable structur, the mutual 
    %                   structure, the classical correlation
    %                   structure), the bilinear odd beads (same order),
    %                   the classical slice of connected beads (even, odd), the
    %                   target qubit, and the difference Q-beads 1 and 2
    %            EntNorm: The bipartite entanglement norm of the state

    % Handle default parameters
    if nargin < 5
        OptStructure = 0;
    end
    if nargin < 6
        MeasurementDirection = [];
    end

    % Check that the bead types are a non-empty list of strings or
    % character arrays
    assert((isstring(BeadTypes) || ischar(BeadTypes) || iscellstr(BeadTypes)) && ~isempty(BeadTypes), ...
           "GetBeads2Qubits:invalidBeadTypes", ...
           "BeadTypes must be a non-empty list of strings or character arrays.");

    % Determine which beads are required for the requested bead types
    BeadTypes = string(BeadTypes);

    % The separable structure (Sigma beads) is required for the mutual, the
    % separable, and all (non)classical beads
    NeedSigma = any(ismember(BeadTypes, ["M", "Sigma", "M_Cl", "Sigma_Cl", "D", "Sigma_D"]));

    % The measurement projectors and the structures after the measurement
    % are required for the classical beads, the nonclassical and the difference
    % Q-beads
    NeedProj = any(ismember(BeadTypes, ["M_Cl", "Sigma_Cl", "D", "Sigma_D"]));

    % The mutual structure (M beads) is required for the mutual beads and
    % for the projector analysis
    NeedM = any(ismember(BeadTypes, ["M", "M_Cl", "Sigma_Cl", "D", "Sigma_D"]));

    % Grid points for the beads, parameters are required for SH_Defs_Real
    res_theta = GridPoints(1);
    res_phi = GridPoints(2);

    % Lisa basis decomposition and mapping for 2 qubits
    % Load definitions of real Lisa tensors and (real and complex) spherical harmonics
    SH_Defs_Real;
    struc1 = load('LiSAReal.mat');
    LiSAReal = struc1.LiSAReal;
    struc2 = load('LiSAScalings.mat');
    LiSAScalings = struc2.LiSAScalings;

    % How many pure basis states contribute to mixed quantum state
    states = size(PureStates,1);
    dim = size(PureStates{1},1);

    % Initialize beads_lists for collecting contributions from all pure states
    % and weigh them by norm squared of pure states

    % All beads, total correlation
    id_bead_mixed = zeros(size(Y00));
    bead1_mixed = zeros(size(Y00));
    bead2_mixed = zeros(size(Y00));
    bilinear_even_mixed_T = zeros(size(Y00));
    bilinear_odd_mixed_T = zeros(size(Y00));

    % Bilinear terms, connected
    bilinear_even_mixed_E = zeros(size(Y00));
    bilinear_odd_mixed_E = zeros(size(Y00));

    % Bilinear terms, compound
    bilinear_even_mixed_C = zeros(size(Y00));
    bilinear_odd_mixed_C = zeros(size(Y00));

    % Beads of the separable structure, classical and nonclassical beads, mutual
    % beads, classical mutual and connected beads, and difference Q-beads,
    % these are only computed if required and are initialized with zeros
    bilinear_even_mixed_Sigma = zeros(size(Y00));
    bilinear_odd_mixed_Sigma = zeros(size(Y00));
    bilinear_even_SigmaCL = zeros(size(Y00));
    bilinear_odd_SigmaCL = zeros(size(Y00));
    bilinear_even_SigmaD = zeros(size(Y00));
    bilinear_odd_SigmaD = zeros(size(Y00));
    bilinear_even_M = zeros(size(Y00));
    bilinear_odd_M = zeros(size(Y00));
    bilinear_even_MCl = zeros(size(Y00));
    bilinear_odd_MCl = zeros(size(Y00));
    bilinear_even_ECl = zeros(size(Y00));
    bilinear_odd_ECl = zeros(size(Y00));
    delta_Q_1 = zeros(size(Y00));
    delta_Q_2 = zeros(size(Y00));

    rho_compound_mixed = zeros(dim,dim);
    rho_connected_ = zeros(dim,dim);

    EntNorm = 0;

    % Loop over all pure states of the decomposition
    for i = 1:states

        % Weight of the pure state and its density operator
        norm_state = norm(PureStates{i});
        weight = norm_state^2;
        rho_pure = (PureStates{i}/norm_state)*(PureStates{i}/norm_state)';

        % First calculate all beads for Total Correlation

        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_T = {};
        [scaled_coeff_T] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_T,rho_pure);

        % Select spherical harmonics to the corresponding coefficients

        % Identity bead
        id_bead = scaled_coeff_T{1}*Y00;
        id_bead_mixed = id_bead_mixed+weight*id_bead;

        % Bead 1 & 2
        bead1 = scaled_coeff_T{2}*Y1n1+scaled_coeff_T{3}*Y10+scaled_coeff_T{4}*Y11;
        bead2 = scaled_coeff_T{5}*Y1n1+scaled_coeff_T{6}*Y10+scaled_coeff_T{7}*Y11;

        bead1_mixed = bead1_mixed+weight*bead1;
        bead2_mixed = bead2_mixed+weight*bead2;

        % Bilinear terms

        % T_00 even and T_2x 5 terms
        bilinear_even_T = scaled_coeff_T{8}*Y00+ scaled_coeff_T{12}*Y2n2+...
                          scaled_coeff_T{13}*Y2n1+scaled_coeff_T{14}*Y20+...
                          scaled_coeff_T{15}*Y21+scaled_coeff_T{16}*Y22;

        bilinear_odd_T = scaled_coeff_T{9}*Y1n1+scaled_coeff_T{10}*Y10+scaled_coeff_T{11}*Y11;

        bilinear_even_mixed_T = bilinear_even_mixed_T+weight*bilinear_even_T;
        bilinear_odd_mixed_T = bilinear_odd_mixed_T+weight*bilinear_odd_T;

        % Calculate rho connected and rho compound for each pure state
        [rho_connected,rho_compound] = rho_connected_pure(rho_pure,dim);

        % Entanglement norm of the pure state, weighted with its probability
        rho_E12 = removeIdLinear(rho_connected);
        EntNorm = EntNorm + weight*getBipEntanglementNorm(4/3*real(trace(rho_E12*rho_E12)));

        % Collect all rho_compounds from the pure states and build the weighted
        % sum for further analysis of the compound part later
        rho_compound_mixed = rho_compound_mixed+weight*rho_compound;
        rho_connected_ = rho_connected_+weight*rho_connected;

        % Calculate Connected Bilinear Terms E

        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_E = {};
        [scaled_coeff_E] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_E,rho_connected);

        % Select spherical harmonics to the corresponding coefficients

        % Bilinear terms

        % T_00 even and T_2x 5 terms
        bilinear_even_E = scaled_coeff_E{8}*Y00+ scaled_coeff_E{12}*Y2n2+...
                          scaled_coeff_E{13}*Y2n1+scaled_coeff_E{14}*Y20+...
                          scaled_coeff_E{15}*Y21+scaled_coeff_E{16}*Y22;

        bilinear_odd_E = scaled_coeff_E{9}*Y1n1+scaled_coeff_E{10}*Y10+scaled_coeff_E{11}*Y11;

        bilinear_even_mixed_E = bilinear_even_mixed_E+weight*bilinear_even_E;
        bilinear_odd_mixed_E = bilinear_odd_mixed_E+weight*bilinear_odd_E;

        % Bilinear terms for rho compound C

        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_C = {};
        [scaled_coeff_C] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_C,rho_compound);

        % Select spherical harmonics to the corresponding coefficients

        % Bilinear terms

        % T_00 even and T_2x 5 terms
        bilinear_even_C = scaled_coeff_C{8}*Y00+ scaled_coeff_C{12}*Y2n2+...
                          scaled_coeff_C{13}*Y2n1+scaled_coeff_C{14}*Y20+...
                          scaled_coeff_C{15}*Y21+scaled_coeff_C{16}*Y22;

        bilinear_odd_C = scaled_coeff_C{9}*Y1n1+scaled_coeff_C{10}*Y10+scaled_coeff_C{11}*Y11;

        bilinear_even_mixed_C = bilinear_even_mixed_C+weight*bilinear_even_C;
        bilinear_odd_mixed_C = bilinear_odd_mixed_C+weight*bilinear_odd_C;
    end

    % Further analysis of rho_compound_mixed

    % Put together original rho for further analysis
    rho_original = zeros(dim,dim);

    for l = 1:states
        rho_original = rho_original+PureStates{l}*PureStates{l}';
    end

    % Beads of the separable structure, only required for the corresponding
    % bead types
    if NeedSigma
        % Get rho after subtracting non correlated part (redundant information)
        rho_sigma = compound_analysis(rho_original,rho_compound_mixed,zeros(dim,dim),dim);

        % Beads of the separable structure rho_sigma

        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_Sigma = {};
        [scaled_coeff_Sigma] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_Sigma,rho_sigma);

        % Select spherical harmonics to the corresponding coefficients

        % Bilinear terms

        % T_00 even and T_2x 5 terms
        bilinear_even_Sigma = scaled_coeff_Sigma{8}*Y00+ scaled_coeff_Sigma{12}*Y2n2+...
                           scaled_coeff_Sigma{13}*Y2n1+scaled_coeff_Sigma{14}*Y20+...
                           scaled_coeff_Sigma{15}*Y21+scaled_coeff_Sigma{16}*Y22;

        bilinear_odd_Sigma = scaled_coeff_Sigma{9}*Y1n1+scaled_coeff_Sigma{10}*Y10+scaled_coeff_Sigma{11}*Y11;

        bilinear_even_mixed_Sigma = bilinear_even_Sigma;
        bilinear_odd_mixed_Sigma = bilinear_odd_Sigma;
    end

    % Mutual structure and its beads, only required for the corresponding
    % bead types
    if NeedM
        rho_M = rho_mutualinformation(rho_original);

        % Beads of the mutual structure
        scaled_coefficients_M = {};
        [scaled_coeff_M] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_M,rho_M);

        % Select spherical harmonics to the corresponding coefficients

        % Bilinear terms

        % T_00 even and T_2x 5 terms
        bilinear_even_M = scaled_coeff_M{8}*Y00+ scaled_coeff_M{12}*Y2n2+...
                          scaled_coeff_M{13}*Y2n1+scaled_coeff_M{14}*Y20+...
                          scaled_coeff_M{15}*Y21+scaled_coeff_M{16}*Y22;

        bilinear_odd_M = scaled_coeff_M{9}*Y1n1+scaled_coeff_M{10}*Y10+scaled_coeff_M{11}*Y11;
    end

    % Further analysis of rho_sigma and rho_connectrf to separate the classical from
    % nonclassical correlations: measurement projectors, classical, nonclassical, and
    % difference Q-beads, only required for the corresponding bead types
    if NeedProj
        % Apply the optimal measurement projectors to the structures
        [rho_sigmaCl,rho_MCl,rho_sigmaD, rho_ECl, drho1, drho2] = find_projectors_bipartite(rho_original,rho_sigma,rho_M,rho_connected_,TargetQubit, OptStructure, MeasurementDirection);

        % CME beads
        scaled_coefficients_ECl = {};
        [scaled_coeff_ECl] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_ECl,rho_ECl);

        % Select spherical harmonics to the corresponding coefficients

        % Bilinear terms

        % T_00 even and T_2x 5 terms
        bilinear_even_ECl = scaled_coeff_ECl{8}*Y00+ scaled_coeff_ECl{12}*Y2n2+...
                            scaled_coeff_ECl{13}*Y2n1+scaled_coeff_ECl{14}*Y20+...
                            scaled_coeff_ECl{15}*Y21+scaled_coeff_ECl{16}*Y22;

        bilinear_odd_ECl = scaled_coeff_ECl{9}*Y1n1+scaled_coeff_ECl{10}*Y10+scaled_coeff_ECl{11}*Y11;

        % Beads of the classical part of the separable structure

        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_SigmaCL = {};
        [scaled_coeff_SigmaCL] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_SigmaCL,rho_sigmaCl);

        % Select spherical harmonics to the corresponding coefficients

        % Bilinear terms

        % T_00 even and T_2x 5 terms
        bilinear_even_SigmaCL = scaled_coeff_SigmaCL{8}*Y00+ scaled_coeff_SigmaCL{12}*Y2n2+...
                           scaled_coeff_SigmaCL{13}*Y2n1+scaled_coeff_SigmaCL{14}*Y20+...
                           scaled_coeff_SigmaCL{15}*Y21+scaled_coeff_SigmaCL{16}*Y22;

        bilinear_odd_SigmaCL = scaled_coeff_SigmaCL{9}*Y1n1+scaled_coeff_SigmaCL{10}*Y10+scaled_coeff_SigmaCL{11}*Y11;

        % Beads of the nonclassical correlation structure

        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_SigmaD = {};
        [scaled_coeff_SigmaD] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_SigmaD,rho_sigmaD);

        % Select spherical harmonics to the corresponding coefficients

        % Bilinear terms

        % T_00 even and T_2x 5 terms
        bilinear_even_SigmaD = scaled_coeff_SigmaD{8}*Y00+ scaled_coeff_SigmaD{12}*Y2n2+...
                          scaled_coeff_SigmaD{13}*Y2n1+scaled_coeff_SigmaD{14}*Y20+...
                          scaled_coeff_SigmaD{15}*Y21+scaled_coeff_SigmaD{16}*Y22;

        bilinear_odd_SigmaD = scaled_coeff_SigmaD{9}*Y1n1+scaled_coeff_SigmaD{10}*Y10+scaled_coeff_SigmaD{11}*Y11;

        % Classical mutual information beads
        scaled_coefficients_MCl = {};
        [scaled_coeff_MCl] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_MCl,rho_MCl);

        % Select spherical harmonics to the corresponding coefficients

        % Bilinear terms

        % T_00 even and T_2x 5 terms
        % Suppression factor of the classical mutual information beads
        % (1: no suppression)
        supp = 1;
        bilinear_even_MCl = supp*scaled_coeff_MCl{8}*Y00+ supp*scaled_coeff_MCl{12}*Y2n2+...
                           supp*scaled_coeff_MCl{13}*Y2n1+supp*scaled_coeff_MCl{14}*Y20+...
                           supp*scaled_coeff_MCl{15}*Y21+supp*scaled_coeff_MCl{16}*Y22;

        bilinear_odd_MCl = supp*scaled_coeff_MCl{9}*Y1n1+supp*scaled_coeff_MCl{10}*Y10+supp*scaled_coeff_MCl{11}*Y11;

        % Difference Q-Beads, the local state of each qubit after the
        % measurement (with a maximally mixed other qubit) minus the original
        % single-qubit beads
        scaled_coefficients_drho1 = {};
        [scaled_coeff_drho1] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_drho1, kron(drho1, eye(2)/2));

        scaled_coefficients_drho2 = {};
        [scaled_coeff_drho2] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_drho2, kron(eye(2)/2, drho2));

        delta_Q_1 = scaled_coeff_drho1{2}*Y1n1+scaled_coeff_drho1{3}*Y10+scaled_coeff_drho1{4}*Y11 - bead1_mixed;
        delta_Q_2 = scaled_coeff_drho2{5}*Y1n1+scaled_coeff_drho2{6}*Y10+scaled_coeff_drho2{7}*Y11 - bead2_mixed;
    end

    % Collect beads for output according to EVEN and then ODD case
    Beads = {id_bead_mixed,bead1_mixed,bead2_mixed,...
             bilinear_even_mixed_T,bilinear_even_mixed_E,...
             bilinear_even_mixed_C,bilinear_even_mixed_Sigma,....
             bilinear_even_SigmaCL,bilinear_even_SigmaD,...
             bilinear_even_M,bilinear_even_MCl,...
             bilinear_odd_mixed_T,bilinear_odd_mixed_E,...
             bilinear_odd_mixed_C,bilinear_odd_mixed_Sigma,...
             bilinear_odd_SigmaCL,bilinear_odd_SigmaD,...
             bilinear_odd_M,bilinear_odd_MCl,...
             bilinear_even_ECl,bilinear_odd_ECl ...
             TargetQubit, delta_Q_1, delta_Q_2};
end
