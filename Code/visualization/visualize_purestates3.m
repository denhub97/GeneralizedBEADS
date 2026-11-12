% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [] = visualize_purestates3(PureStates, CameraAngle, BeadType, Id, Lines, Weights, Tangles, EntArcs, GridPoints, Recursive, Path, FontSize)

    %!!! visualize_purestates3 visualizes the BEADS representations of all
    %    pure states constituting a decomposition of a mixed 3-qubit state,
    %    each pure state in an individual subplot (up to 8 per figure, further
    %    states are plotted into additional figures), for one bead type (T,
    %    E, or C). The coordinate axes are drawn in the first subplot.
    %    Optional separator lines, the weights of the pure states, and the
    %    two-tangles and three-tangle can be added
    %
    %    inputs: PureStates: cell array of the contributing pure states (the
    %                        weights are contained in the state vectors,
    %                        i.e., sqrt(p)*|state>)
    %            CameraAngle: the camera angle [rotation, elevation]
    %            BeadType: the type of correlation function beads ("T", "E",
    %                      "C")
    %            Id: specify whether id bead shall be drawn (boolean)
    %            Lines: determines whether separator lines are drawn between
    %                   the individual subplots (boolean)
    %            Weights: determines whether the weight of each pure state
    %                     is shown in the top-left corner of each subplot
    %                     (boolean)
    %            Tangles: determines whether the two-tangles and the
    %                     three-tangle of each pure state are shown in the
    %                     top-left corner (boolean)
    %            EntArcs: draw connector bars scaled by the specific
    %                     tangles (boolean)
    %            GridPoints: dimensions of meshgrid for plotting [theta, phi]
    %            Recursive: determines whether the reduced marginals shall
    %                       be optimized individually for plotting (boolean)
    %            Path: path to where the final image shall be saved, the
    %                  figures are numbered if there are more than 8 states
    %            FontSize: font size for all labels and annotations
    %
    %    outputs: none, the figures are saved to Path

    % Grid points for the beads
    res_theta = GridPoints(1);
    res_phi = GridPoints(2);

    % Load all definitions and external data (LiSA operators, spherical harm.)
    SH_Defs_Real;
    struc1 = load('LiSAReal.mat');
    LiSAReal = struc1.LiSAReal;
    struc2 = load('LiSAScalings.mat');
    LiSAScalings = struc2.LiSAScalings;

    % How many pure basis states contribute to the original mixed quantum state
    states = size(PureStates,1);
    dim = size(PureStates{1},1);

    % Layout parameters
    label_shift = 0.7;
    label_fontsize = FontSize;

    tolerance = 1e-2;

    % The basic beads positions
    beads_positions = [-7, -4.5, 0;  % Id
                       -6,  0,   0;  % Q1
                       0,  0,   0;  % Q2
                       6,  0,   0;  % Q3
                       -3,  2,   0;  % {1,2}_even
                       -3,  4.5, 0;  % {1,2}_odd
                       0, -3,   0;  % {1,3}_even
                       0, -5.5, 0;  % {1,3}_odd
                       3,  2,   0;  % {2,3}_even
                       3,  4.5, 0;  % {2,3}_odd
                       0,  7.5, 0;  % Tau1_odd
                       -5, 10.5, 0;  % Tau2_even
                       -5, 13.5, 0;  % Tau2_odd
                       0, 10.5, 0;  % Tau3_even
                       0, 13.5, 0;  % Tau3_odd
                       5, 10.5, 0]; % Tau4_even

    % Determine the viewpoint
    % Top view
    if CameraAngle(1) == 0 && CameraAngle(2) == 90
        viewpoint = "top";

        % Label positions
        label_positions = [label_shift-7,        label_shift-4.5,       0; % Id
                           label_shift*1.5-6,    0,                     0; % Q1
                           label_shift*1.5,      0,                     0; % Q2
                           label_shift*1.5+6,    0,                     0; % Q3
                           label_shift-3,        label_shift+2,         0; % {1,2}_even
                           label_shift-3,        label_shift+4.5,       0; % {1,2}_odd
                           label_shift,          label_shift-3,         0; % {1,3}_even
                           label_shift,          label_shift-5.5,       0; % {1,3}_odd
                           label_shift+3,        label_shift+2,         0; % {2,3}_even
                           label_shift+3,        label_shift+4.5,       0; % {2,3}_odd
                           label_shift,          label_shift+7.5,       0; % Tau1_odd
                           label_shift-5,        label_shift*1.5+10.5,  0; % Tau2_even
                           label_shift-5,        label_shift*1.5+13.5,  0; % Tau2_odd
                           label_shift,          label_shift*1.5+10.5,  0; % Tau3_even
                           label_shift,          label_shift*1.5+13.5,  0; % Tau3_odd
                           label_shift+5,        label_shift*1.5+10.5,  0];% Tau4_even

    % Oblique views
    elseif CameraAngle(1) == 54 && CameraAngle(2) == 45
        viewpoint = "oblique";

        if EntArcs
            % Label positions for the case that connector lines are enabled

            label_positions = [label_shift*1.2-7,       label_shift*1.2-4.5,       0; % Id
                               label_shift*2-6,         0,                        -2.2*label_shift; % Q1
                               label_shift*2,           0,                        -2.2*label_shift; % Q2
                               label_shift*2+6,         0,                        -2.2*label_shift; % Q3
                               label_shift*3-10.5,      label_shift*(-0.2)+1,     -0.4; % {1,2}_even
                               label_shift*(-1.7)-6.5,  4.5,                        -0.5; % {1,2}_odd
                               label_shift*(-6),        label_shift-8,           -0.5; % {1,3}_even
                               label_shift*2.3,         -5.5,                      -2.2*label_shift; % {1,3}_odd
                               label_shift*2.5+4.5,     label_shift*0.2+2.5,      -2*label_shift; % {2,3}_even
                               label_shift+4,           label_shift*1.2+4.5,      -2.2*label_shift; % {2,3}_odd
                               label_shift*3.5-5.5,     label_shift*0.6+6,        -0.2; % Tau1_odd
                               label_shift*(-1)-7.8,    label_shift*(-1.5)+8,      0.3; % Tau2_even
                               label_shift-5,           label_shift*1.5+13.5,      0; % Tau2_odd
                               label_shift*2,           11,                       -1.5; % Tau3_even
                               label_shift,             label_shift*1.5+13.5,      0; % Tau3_odd
                               label_shift+4.5,         label_shift*1.5+10.5,     -3];% Tau4_even

        else
            % Label positions
            label_positions = [label_shift*1.2-7,       label_shift*1.2-4.5,       0; % Id
                               label_shift*2-6,         0,                        -2.2*label_shift; % Q1
                               label_shift*2,           0,                        -2.2*label_shift; % Q2
                               label_shift*2+6,         0,                        -2.2*label_shift; % Q3
                               label_shift*3-3,         label_shift*(-0.2)+2,     -0.4; % {1,2}_even
                               label_shift*(-1.7)-5.7,  1,                         0; % {1,2}_odd
                               label_shift*(-6.0),      label_shift-7.5,           0; % {1,3}_even
                               label_shift*2.3,        -5.5,                      -2.2*label_shift; % {1,3}_odd
                               label_shift*2.5+3,       label_shift*0.2+2,        -2*label_shift; % {2,3}_even
                               label_shift+3.5,         label_shift*1.2+4.5,      -2.2*label_shift; % {2,3}_odd
                               label_shift*2.1-5,       label_shift*0.6+5.5,      -0.4; % Tau1_odd
                               label_shift*(-1)-7.8,    label_shift*(-1.5)+8,      0; % Tau2_even
                               label_shift-5,           label_shift*1.5+13.5,      0; % Tau2_odd
                               label_shift*2,           11,                       -1.5; % Tau3_even
                               label_shift,             label_shift*1.5+13.5,      0; % Tau3_odd
                               label_shift+4.5,         label_shift*1.5+10.5,     -3];% Tau4_even
        end

    % Other views we don't have labels
    else
        viewpoint = "none";
    end

    % Determine the beads labels
    switch BeadType
        % Total correlation function
        case "T"
            labels = {'Id','$Q_1$','$Q_2$','$Q_3$','$T\{1,2\}_{even}$',...
                      '$T\{1,2\}_{odd}$','$T\{1,3\}_{even}$','$T\{1,3\}_{odd}$',...
                      '$T\{2,3\}_{even}$','$T\{2,3\}_{odd}$',...
                      '$T \{\tau_1\}_{odd}$','$T \{\tau_2\}_{even}$',...
                      '$T \{\tau_2\}_{odd}$','$T \{\tau_3\}_{even}$',...
                      '$T \{\tau_3\}_{odd}$','$T \{\tau_4\}_{even}$'};

        % Connected (entanglement-based) correlation function
        case "E"
            labels = {'Id','$Q_1$','$Q_2$','$Q_3$','$E\{1,2\}_{even}$',...
                      '$E\{1,2\}_{odd}$','$E\{1,3\}_{even}$','$E\{1,3\}_{odd}$',...
                      '$E\{2,3\}_{even}$','$E\{2,3\}_{odd}$',...
                      '$E \{\tau_1\}_{odd}$','$E \{\tau_2\}_{even}$',...
                      '$E \{\tau_2\}_{odd}$','$E \{\tau_3\}_{even}$',...
                      '$E \{\tau_3\}_{odd}$','$E \{\tau_4\}_{even}$'};

        % Compound correlation function
        case "C"
            labels = {'Id','$Q_1$','$Q_2$','$Q_3$','$C\{1,2\}_{even}$',...
                      '$C\{1,2\}_{odd}$','$C\{1,3\}_{even}$','$C\{1,3\}_{odd}$',...
                      '$C\{2,3\}_{even}$','$C\{2,3\}_{odd}$',...
                      '$C \{\tau_1\}_{odd}$','$C \{\tau_2\}_{even}$',...
                      '$C \{\tau_2\}_{odd}$','$C \{\tau_3\}_{even}$',...
                      '$C \{\tau_3\}_{odd}$','$C \{\tau_4\}_{even}$'};
    end

    % Coordinate axes layout parameters
    cylinder_radius = 0.05;
    cone_height = 0.2;
    cone_radius = 0.1;
    axis_position = [-9.5, 0, 0];
    if viewpoint == "oblique" && EntArcs
        axis_position = [-11, 0, 0];
    end

    % Permutations to get all 3D coordinate axes
    axis_perm = {[1 2 3], [1 3 2], [3 2 1]};
    axis_colors = [0  0  0;
                   .9 .8  0;
                   1  0  0];
    axis_labels = {'z', 'y', 'x'};

    % Define spherical coordinate mesh for plotting
    theta = linspace(0,pi,GridPoints(1));
    phi = linspace(0,2*pi,GridPoints(2));
    [Phi,Theta] = meshgrid(phi,theta);
    radius=1;

    % Cartesian from spherical coordinates for 3D meshgrid
    X = cos(Phi).*sin(Theta).*radius;
    Y = sin(Phi).*sin(Theta).*radius;
    Z = cos(Theta)*radius;

    % The pure states are plotted in groups of 8 states per figure, count is
    % the number of states which have been plotted already
    count = 0;
    while count < states
        path_ = regexprep(Path, '\.(png|jpg|pdf|tiff)$', '_pure_decomp' + BeadType + num2str(count) + '.$1');

        % Create the figure and subplot
        fig = figure;
        if isscalar(PureStates)
            nCols = 1;
        else
            nCols = 2;
        end
        nRows = min(4,ceil(states/nCols));
        fig.Position = [0,0,700,600*nRows/nCols];
        gap = [0.05 0.];
        marg_h = [.03 .02];
        marg_w = [.0 .0];
        ha = tight_subplot(nRows, nCols, gap, marg_h, marg_w);

        % Zero array, required for color assignment
        blank = 0*X;

        % List of probabilities, initialize empty, filled in next step
        probabilities = zeros(states,1);
        two_tangles = zeros(states,3);
        three_tangles = zeros(states,1);

        ent_norm_123 = zeros(numel(states),1);

        % Plot each pure state in an individual subplot
        maxidx = min(8, states-count);
        for i = 1:maxidx
            % Norm of the weighted pure state
            norm_state = norm(PureStates{count+i});
            probabilities(i) = norm_state^2;

            % Renormalize pure state and calculate density operator
            rho_pure = (PureStates{count+i}/norm_state)*(PureStates{count+i}/norm_state)';

            % Get reduced density matrices for the 2 qubits subsystems
            rho_12 = partial_trace(rho_pure,3);
            rho_13 = partial_trace(rho_pure,2);
            rho_23 = partial_trace(rho_pure,1);

            % For the recursive treatment, optimize the decompositions of the two-qubit
            % marginals of the pure state individually
            if Recursive
                [~,states_12] = opt_ent_formation(rho_12);
                [~,states_13] = opt_ent_formation(rho_13);
                [~,states_23] = opt_ent_formation(rho_23);
            end

            % First calculate 2-qubit pairwise concurrence
            C_12 = get_concurrence(rho_12);
            C_13 = get_concurrence(rho_13);
            C_23 = get_concurrence(rho_23);

            % Two tangles are the squared concurrences
            two_tangles(i,:) = [C_12^2, C_13^2, C_23^2];

            % Calculate three tangles
            [~,~,~,~,tau3] = get_three_tangle(rho_pure);
            three_tangles(i) = abs(tau3);

            % First calculate all beads for total correlation function
            % Get scaled coefficients of Lisa tensor basis operators
            scaled_coefficients = {};
            [scaled_coeff]=traverse_and_multiply(LiSAReal{3},LiSAScalings{3},...
                scaled_coefficients,rho_pure);

            % Calculate rho connected and rho compound for each pure state
            [rho_connected,rho_compound] = rho_connected_pure(rho_pure,dim);

            % Calculate the tripartite entanglement norm
            [~,~,~,rho_E123] = split_op_tripartite(rho_connected);
            rho_E123 = removeIdLinear(rho_E123);
            % Normalize to 1 (GHZ state has maximum value 1 accordingly)
            ent_norm_123(i) = getEntanglementNorm(2*real(trace(rho_E123*rho_E123)));

            % Connected (entanglement-based) correlation function
            % Get scaled coefficients of Lisa tensor basis
            scaled_coefficients_E = {};
            [scaled_coeff_E]=traverse_and_multiply(LiSAReal{3},LiSAScalings{3},...
                scaled_coefficients_E,rho_connected);

            % Compound correlation function
            % Get scaled coefficients of Lisa tensor basis
            scaled_coefficients_C = {};
            [scaled_coeff_C] = traverse_and_multiply(LiSAReal{3},LiSAScalings{3},...
                scaled_coefficients_C,rho_compound);

            % Scale spherical harmonics by corresponding coefficients
            % Id bead
            id_bead = scaled_coeff{1}*Y00;

            % Q-Beads 1 to 3
            bead1 = scaled_coeff{2}*Y1n1 + scaled_coeff{3}*Y10 + scaled_coeff{4}*Y11;
            bead2 = scaled_coeff{5}*Y1n1 + scaled_coeff{6}*Y10 + scaled_coeff{7}*Y11;
            bead3 = scaled_coeff{8}*Y1n1 + scaled_coeff{9}*Y10 + scaled_coeff{10}*Y11;

            % Trilinear terms, total correlation function
            % Scale spherical harmonics by corresponding coefficients
            % j = 1, m = [-1, ..., 1] and j = 3, m = [-3, ..., 3] terms
            trilin_tau1_odd = scaled_coeff{38}*Y1n1 + scaled_coeff{39}*Y10 +...
                scaled_coeff{40}*Y11 + scaled_coeff{41}*Y3n3 + scaled_coeff{42}*Y3n2 +...
                scaled_coeff{43}*Y3n1 + scaled_coeff{44}*Y30 + ...
                scaled_coeff{45}*Y31 + scaled_coeff{46}*Y32 + scaled_coeff{47}*Y33;

            % j = 2, m = [-2, ..., 2] terms
            trilin_tau2_even = scaled_coeff{51}*Y2n2 + scaled_coeff{52}*Y2n1 + scaled_coeff{53}*Y20 +...
                scaled_coeff{54}*Y21 + scaled_coeff{55}*Y22;

            % j = 1, m = [-1, ..., 1] terms
            trilin_tau2_odd = scaled_coeff{48}*Y1n1 + scaled_coeff{49}*Y10 + scaled_coeff{50}*Y11;

            % j = 1, m = [-1, ..., 1] terms
            trilin_tau3_odd = scaled_coeff{56}*Y1n1 + scaled_coeff{57}*Y10 + scaled_coeff{58}*Y11;
            % j = 2, m = [-2, ..., 2] terms
            trilin_tau3_even = scaled_coeff{59}*Y2n2 + scaled_coeff{60}*Y2n1 + scaled_coeff{61}*Y20 +...
                scaled_coeff{62}*Y21 + scaled_coeff{63}*Y22;
            % j = 0, m = 0
            trilin_tau4_even = scaled_coeff{64}*Y00;

            % Trilinear terms, connected (entanglement-based) correlation function
            % Scale spherical harmonics by corresponding coefficients
            % j = 1, m = [-1, ..., 1] and j = 3, m = [-3, ..., 3] terms
            trilin_tau1_odd_E = scaled_coeff_E{38}*Y1n1 + scaled_coeff_E{39}*Y10 +...
                scaled_coeff_E{40}*Y11 + scaled_coeff_E{41}*Y3n3 + scaled_coeff_E{42}*Y3n2 +...
                scaled_coeff_E{43}*Y3n1 + scaled_coeff_E{44}*Y30 + ...
                scaled_coeff_E{45}*Y31 + scaled_coeff_E{46}*Y32 + scaled_coeff_E{47}*Y33;

            % j = 2, m = [-2, ..., 2] terms
            trilin_tau2_even_E = scaled_coeff_E{51}*Y2n2 + scaled_coeff_E{52}*Y2n1 + scaled_coeff_E{53}*Y20 +...
                scaled_coeff_E{54}*Y21 + scaled_coeff_E{55}*Y22;

            % j = 1, m = [-1, ..., 1] terms
            trilin_tau2_odd_E = scaled_coeff_E{48}*Y1n1 + scaled_coeff_E{49}*Y10 + scaled_coeff_E{50}*Y11;

            % j = 1, m = [-1, ..., 1] terms
            trilin_tau3_odd_E = scaled_coeff_E{56}*Y1n1 + scaled_coeff_E{57}*Y10 + scaled_coeff_E{58}*Y11;
            % j = 2, m = [-2, ..., 2] terms
            trilin_tau3_even_E = scaled_coeff_E{59}*Y2n2 + scaled_coeff_E{60}*Y2n1 + scaled_coeff_E{61}*Y20 +...
                scaled_coeff_E{62}*Y21 + scaled_coeff_E{63}*Y22;
            % j = 0, m = 0
            trilin_tau4_even_E = scaled_coeff_E{64}*Y00;

            % Trilinear terms, compound correlation function
            % Scale spherical harmonics by corresponding coefficients
            % j = 1, m = [-1, ..., 1] and j = 3, m = [-3, ..., 3] terms
            trilin_tau1_odd_C = scaled_coeff_C{38}*Y1n1 + scaled_coeff_C{39}*Y10 +...
                scaled_coeff_C{40}*Y11 + scaled_coeff_C{41}*Y3n3 + scaled_coeff_C{42}*Y3n2 +...
                scaled_coeff_C{43}*Y3n1 + scaled_coeff_C{44}*Y30 + ...
                scaled_coeff_C{45}*Y31 + scaled_coeff_C{46}*Y32 + scaled_coeff_C{47}*Y33;

            % j = 2, m = [-2, ..., 2] terms
            trilin_tau2_even_C = scaled_coeff_C{51}*Y2n2 + scaled_coeff_C{52}*Y2n1 + scaled_coeff_C{53}*Y20 +...
                scaled_coeff_C{54}*Y21 + scaled_coeff_C{55}*Y22;

            % j = 1, m = [-1, ..., 1] terms
            trilin_tau2_odd_C = scaled_coeff_C{48}*Y1n1 + scaled_coeff_C{49}*Y10 + scaled_coeff_C{50}*Y11;

            % j = 1, m = [-1, ..., 1] terms
            trilin_tau3_odd_C = scaled_coeff_C{56}*Y1n1 + scaled_coeff_C{57}*Y10 + scaled_coeff_C{58}*Y11;
            % j = 2, m = [-2, ..., 2] terms
            trilin_tau3_even_C = scaled_coeff_C{59}*Y2n2 + scaled_coeff_C{60}*Y2n1 + scaled_coeff_C{61}*Y20 +...
                scaled_coeff_C{62}*Y21 + scaled_coeff_C{63}*Y22;
            % j = 0, m = 0
            trilin_tau4_even_C = scaled_coeff_C{64}*Y00;

            % Bilinear beads of the marginals, either directly from the pure state
            % (not recursive) or from the optimized decompositions of the marginals
            % (recursive)
            if Recursive == false
                % Bilinear terms, total correlation function
                % Scale spherical harmonics by corresponding coefficients
                % j = 0, m = 0 and j = 2, m = [-2, ..., 2] terms
                bilin12_even = scaled_coeff{11}*Y00 + scaled_coeff{15}*Y2n2 +...
                    scaled_coeff{16}*Y2n1 + scaled_coeff{17}*Y20 +...
                    scaled_coeff{18}*Y21 + scaled_coeff{19}*Y22;
                bilin13_even = scaled_coeff{20}*Y00 + scaled_coeff{24}*Y2n2 +...
                    scaled_coeff{25}*Y2n1 + scaled_coeff{26}*Y20 +...
                    scaled_coeff{27}*Y21 + scaled_coeff{28}*Y22;
                bilin23_even = scaled_coeff{29}*Y00 + scaled_coeff{33}*Y2n2 +...
                    scaled_coeff{34}*Y2n1 + scaled_coeff{35}*Y20 +...
                    scaled_coeff{36}*Y21 + scaled_coeff{37}*Y22;

                % j = 1, m = [-1, ..., 1] terms
                bilin12_odd = scaled_coeff{12}*Y1n1 + scaled_coeff{13}*Y10 + scaled_coeff{14}*Y11;
                bilin13_odd = scaled_coeff{21}*Y1n1 + scaled_coeff{22}*Y10 + scaled_coeff{23}*Y11;
                bilin23_odd = scaled_coeff{30}*Y1n1 + scaled_coeff{31}*Y10 + scaled_coeff{32}*Y11;

                % Bilinear terms, connected (entanglement-based) correlation function
                % Scale spherical harmonics by corresponding coefficients
                % j = 0, m = 0 and j = 2, m = [-2, ..., 2] terms
                bilin12_even_E = scaled_coeff_E{11}*Y00 + scaled_coeff_E{15}*Y2n2 +...
                    scaled_coeff_E{16}*Y2n1 + scaled_coeff_E{17}*Y20 +...
                    scaled_coeff_E{18}*Y21 + scaled_coeff_E{19}*Y22;
                bilin13_even_E = scaled_coeff_E{20}*Y00 + scaled_coeff_E{24}*Y2n2 +...
                    scaled_coeff_E{25}*Y2n1 + scaled_coeff_E{26}*Y20 +...
                    scaled_coeff_E{27}*Y21 + scaled_coeff_E{28}*Y22;
                bilin23_even_E = scaled_coeff_E{29}*Y00 + scaled_coeff_E{33}*Y2n2 +...
                    scaled_coeff_E{34}*Y2n1 + scaled_coeff_E{35}*Y20 +...
                    scaled_coeff_E{36}*Y21 + scaled_coeff_E{37}*Y22;

                % j = 1, m = [-1, ..., 1] terms
                bilin12_odd_E = scaled_coeff_E{12}*Y1n1 + scaled_coeff_E{13}*Y10 + scaled_coeff_E{14}*Y11;
                bilin13_odd_E = scaled_coeff_E{21}*Y1n1 + scaled_coeff_E{22}*Y10 + scaled_coeff_E{23}*Y11;
                bilin23_odd_E = scaled_coeff_E{30}*Y1n1 + scaled_coeff_E{31}*Y10 + scaled_coeff_E{32}*Y11;

                % Bilinear terms, compound correlation function
                % Scale spherical harmonics by corresponding coefficients
                % j = 0, m = 0 and j = 2, m = [-2, ..., 2] terms
                bilin12_even_C = scaled_coeff_C{11}*Y00 + scaled_coeff_C{15}*Y2n2 +...
                    scaled_coeff_C{16}*Y2n1 + scaled_coeff_C{17}*Y20 +...
                    scaled_coeff_C{18}*Y21 + scaled_coeff_C{19}*Y22;
                bilin13_even_C = scaled_coeff_C{20}*Y00 + scaled_coeff_C{24}*Y2n2 +...
                    scaled_coeff_C{25}*Y2n1 + scaled_coeff_C{26}*Y20 +...
                    scaled_coeff_C{27}*Y21 + scaled_coeff_C{28}*Y22;
                bilin23_even_C = scaled_coeff_C{29}*Y00 + scaled_coeff_C{33}*Y2n2 +...
                    scaled_coeff_C{34}*Y2n1 + scaled_coeff_C{35}*Y20 +...
                    scaled_coeff_C{36}*Y21 + scaled_coeff_C{37}*Y22;

                % j = 1, m = [-1, ..., 1] terms
                bilin12_odd_C = scaled_coeff_C{12}*Y1n1 + scaled_coeff_C{13}*Y10 + scaled_coeff_C{14}*Y11;
                bilin13_odd_C = scaled_coeff_C{21}*Y1n1 + scaled_coeff_C{22}*Y10 + scaled_coeff_C{23}*Y11;
                bilin23_odd_C = scaled_coeff_C{30}*Y1n1 + scaled_coeff_C{31}*Y10 + scaled_coeff_C{32}*Y11;
            else
                bilin12_even = zeros(GridPoints(1), GridPoints(2));
                bilin12_odd = zeros(GridPoints(1), GridPoints(2));
                bilin12_even_E = zeros(GridPoints(1), GridPoints(2));
                bilin12_odd_E = zeros(GridPoints(1), GridPoints(2));
                bilin12_even_C = zeros(GridPoints(1), GridPoints(2));
                bilin12_odd_C = zeros(GridPoints(1), GridPoints(2));
                rho_connected_mixed_12 = zeros(4);
                rho_compound_mixed_12 = zeros(4);

                bilin13_even = zeros(GridPoints(1), GridPoints(2));
                bilin13_odd = zeros(GridPoints(1), GridPoints(2));
                bilin13_even_E = zeros(GridPoints(1), GridPoints(2));
                bilin13_odd_E = zeros(GridPoints(1), GridPoints(2));
                bilin13_even_C = zeros(GridPoints(1), GridPoints(2));
                bilin13_odd_C = zeros(GridPoints(1), GridPoints(2));
                rho_connected_mixed_13 = zeros(4);
                rho_compound_mixed_13 = zeros(4);

                bilin23_even = zeros(GridPoints(1), GridPoints(2));
                bilin23_odd = zeros(GridPoints(1), GridPoints(2));
                bilin23_even_E = zeros(GridPoints(1), GridPoints(2));
                bilin23_odd_E = zeros(GridPoints(1), GridPoints(2));
                bilin23_even_C = zeros(GridPoints(1), GridPoints(2));
                bilin23_odd_C = zeros(GridPoints(1), GridPoints(2));
                rho_connected_mixed_23 = zeros(4);
                rho_compound_mixed_23 = zeros(4);

                for ii = 1:numel(states_12)
                    norm_state = norm(states_12{ii});
                    weight = norm_state^2;
                    rho_pure_12 = (states_12{ii}/norm_state)*(states_12{ii}/norm_state)';

                    % First calculate all bipartite beads for Total Correlation function
                    % Get scaled coefficients of Lisa tensor basis
                    scaled_coefficients = {};
                    [scaled_coeff] = traverse_and_multiply(LiSAReal{2}, LiSAScalings{2}, scaled_coefficients, rho_pure_12);

                    % Bilinear terms
                    % T_00 even and T_2x 5 terms
                    bilin12_even_ = scaled_coeff{8}*Y00 + scaled_coeff{12}*Y2n2+...
                        scaled_coeff{13}*Y2n1+scaled_coeff{14}*Y20+...
                        scaled_coeff{15}*Y21+scaled_coeff{16}*Y22;

                    bilin12_odd_ = scaled_coeff{9}*Y1n1+scaled_coeff{10}*Y10+scaled_coeff{11}*Y11;

                    bilin12_even = bilin12_even + weight*bilin12_even_;
                    bilin12_odd = bilin12_odd + weight*bilin12_odd_;

                    % Calculate rho connected and rho compound for each pure state
                    [rho_connected_12, rho_compound_12] = rho_connected_pure(rho_pure_12, 4);

                    rho_compound_mixed_12 = rho_compound_mixed_12 + rho_compound_12*weight;
                    rho_connected_mixed_12 = rho_connected_mixed_12 + rho_connected_12*weight;

                    % Bilin 12 Connected
                    scaled_coefficients_E = {};
                    [scaled_coeff_E] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_E,rho_connected_12);

                    bilin12_even_E_ = scaled_coeff_E{8}*Y00+ scaled_coeff_E{12}*Y2n2+...
                        scaled_coeff_E{13}*Y2n1+scaled_coeff_E{14}*Y20+...
                        scaled_coeff_E{15}*Y21+scaled_coeff_E{16}*Y22;

                    bilin12_odd_E_ = scaled_coeff_E{9}*Y1n1+scaled_coeff_E{10}*Y10+scaled_coeff_E{11}*Y11;

                    bilin12_even_E = bilin12_even_E + weight*bilin12_even_E_;
                    bilin12_odd_E = bilin12_odd_E + weight*bilin12_odd_E_;

                    % Bilin 12 Compound
                    scaled_coefficients_C = {};
                    [scaled_coeff_C] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_C,rho_compound_12);

                    bilin12_even_C_ = scaled_coeff_C{8}*Y00+ scaled_coeff_C{12}*Y2n2+...
                        scaled_coeff_C{13}*Y2n1+scaled_coeff_C{14}*Y20+...
                        scaled_coeff_C{15}*Y21+scaled_coeff_C{16}*Y22;

                    bilin12_odd_C_ = scaled_coeff_C{9}*Y1n1+scaled_coeff_C{10}*Y10+scaled_coeff_C{11}*Y11;

                    bilin12_even_C = bilin12_even_C+weight*bilin12_even_C_;
                    bilin12_odd_C = bilin12_odd_C+weight*bilin12_odd_C_;
                end

                for ii = 1:numel(states_13)
                    norm_state = norm(states_13{ii});
                    weight = norm_state^2;
                    rho_pure_13 = (states_13{ii}/norm_state)*(states_13{ii}/norm_state)';

                    % First calculate all bipartite beads for Total Correlation function
                    % Get scaled coefficients of Lisa tensor basis
                    scaled_coefficients = {};
                    [scaled_coeff] = traverse_and_multiply(LiSAReal{2}, LiSAScalings{2}, scaled_coefficients, rho_pure_13);

                    % Bilinear terms
                    % T_00 even and T_2x 5 terms
                    bilin13_even_ = scaled_coeff{8}*Y00 + scaled_coeff{12}*Y2n2+...
                        scaled_coeff{13}*Y2n1+scaled_coeff{14}*Y20+...
                        scaled_coeff{15}*Y21+scaled_coeff{16}*Y22;

                    bilin13_odd_ = scaled_coeff{9}*Y1n1+scaled_coeff{10}*Y10+scaled_coeff{11}*Y11;

                    bilin13_even = bilin13_even + weight*bilin13_even_;
                    bilin13_odd = bilin13_odd + weight*bilin13_odd_;

                    % Calculate rho connected and rho compound for each pure state
                    [rho_connected_13, rho_compound_13] = rho_connected_pure(rho_pure_13, 4);

                    rho_compound_mixed_13 = rho_compound_mixed_13 + rho_compound_13*weight;
                    rho_connected_mixed_13 = rho_connected_mixed_13 + rho_connected_13*weight;

                    % Bilin 13 Connected
                    scaled_coefficients_E = {};
                    [scaled_coeff_E] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_E,rho_connected_13);

                    bilin13_even_E_ = scaled_coeff_E{8}*Y00+ scaled_coeff_E{12}*Y2n2+...
                        scaled_coeff_E{13}*Y2n1+scaled_coeff_E{14}*Y20+...
                        scaled_coeff_E{15}*Y21+scaled_coeff_E{16}*Y22;

                    bilin13_odd_E_ = scaled_coeff_E{9}*Y1n1+scaled_coeff_E{10}*Y10+scaled_coeff_E{11}*Y11;

                    bilin13_even_E = bilin13_even_E + weight*bilin13_even_E_;
                    bilin13_odd_E = bilin13_odd_E + weight*bilin13_odd_E_;

                    % Bilin 13 Compound
                    scaled_coefficients_C = {};
                    [scaled_coeff_C] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_C,rho_compound_13);

                    bilin13_even_C_ = scaled_coeff_C{8}*Y00+ scaled_coeff_C{12}*Y2n2+...
                        scaled_coeff_C{13}*Y2n1+scaled_coeff_C{14}*Y20+...
                        scaled_coeff_C{15}*Y21+scaled_coeff_C{16}*Y22;

                    bilin13_odd_C_ = scaled_coeff_C{9}*Y1n1+scaled_coeff_C{10}*Y10+scaled_coeff_C{11}*Y11;

                    bilin13_even_C = bilin13_even_C+weight*bilin13_even_C_;
                    bilin13_odd_C = bilin13_odd_C+weight*bilin13_odd_C_;
                end

                for ii = 1:numel(states_23)
                    norm_state = norm(states_23{ii});
                    weight = norm_state^2;
                    rho_pure_23 = (states_23{ii}/norm_state)*(states_23{ii}/norm_state)';

                    % First calculate all bipartite beads for Total Correlation function
                    % Get scaled coefficients of Lisa tensor basis
                    scaled_coefficients = {};
                    [scaled_coeff] = traverse_and_multiply(LiSAReal{2}, LiSAScalings{2}, scaled_coefficients, rho_pure_23);

                    % Bilinear terms
                    % T_00 even and T_2x 5 terms
                    bilin23_even_ = scaled_coeff{8}*Y00 + scaled_coeff{12}*Y2n2+...
                        scaled_coeff{13}*Y2n1+scaled_coeff{14}*Y20+...
                        scaled_coeff{15}*Y21+scaled_coeff{16}*Y22;

                    bilin23_odd_ = scaled_coeff{9}*Y1n1+scaled_coeff{10}*Y10+scaled_coeff{11}*Y11;

                    bilin23_even = bilin23_even + weight*bilin23_even_;
                    bilin23_odd = bilin23_odd + weight*bilin23_odd_;

                    % Calculate rho connected and rho compound for each pure state
                    [rho_connected_23, rho_compound_23] = rho_connected_pure(rho_pure_23, 4);

                    rho_compound_mixed_23 = rho_compound_mixed_23 + rho_compound_23*weight;
                    rho_connected_mixed_23 = rho_connected_mixed_23 + rho_connected_23*weight;

                    % Bilin 13 Connected
                    scaled_coefficients_E = {};
                    [scaled_coeff_E] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_E,rho_connected_23);

                    bilin23_even_E_ = scaled_coeff_E{8}*Y00+ scaled_coeff_E{12}*Y2n2+...
                        scaled_coeff_E{13}*Y2n1+scaled_coeff_E{14}*Y20+...
                        scaled_coeff_E{15}*Y21+scaled_coeff_E{16}*Y22;

                    bilin23_odd_E_ = scaled_coeff_E{9}*Y1n1+scaled_coeff_E{10}*Y10+scaled_coeff_E{11}*Y11;

                    bilin23_even_E = bilin23_even_E + weight*bilin23_even_E_;
                    bilin23_odd_E = bilin23_odd_E + weight*bilin23_odd_E_;

                    % Bilin 13 Compound
                    scaled_coefficients_C = {};
                    [scaled_coeff_C] = traverse_and_multiply(LiSAReal{2},LiSAScalings{2},scaled_coefficients_C,rho_compound_23);

                    bilin23_even_C_ = scaled_coeff_C{8}*Y00+ scaled_coeff_C{12}*Y2n2+...
                        scaled_coeff_C{13}*Y2n1+scaled_coeff_C{14}*Y20+...
                        scaled_coeff_C{15}*Y21+scaled_coeff_C{16}*Y22;

                    bilin23_odd_C_ = scaled_coeff_C{9}*Y1n1+scaled_coeff_C{10}*Y10+scaled_coeff_C{11}*Y11;

                    bilin23_even_C = bilin23_even_C+weight*bilin23_even_C_;
                    bilin23_odd_C = bilin23_odd_C+weight*bilin23_odd_C_;
                end
            end

            % Determine the beads' spherical function values (beads) and values
            % used for determination of color (beads_E/C)
            switch BeadType
                % Total correlation function
                case "T"
                    beads = {id_bead, bead1, bead2, bead3, ...
                             bilin12_even, bilin12_odd, ...
                             bilin13_even, bilin13_odd, ...
                             bilin23_even, bilin23_odd, ...
                             trilin_tau1_odd, trilin_tau2_even, ...
                             trilin_tau2_odd, trilin_tau3_even, ...
                             trilin_tau3_odd, trilin_tau4_even};

                    beads_E = {blank, blank, blank, blank, ...
                               bilin12_even_E, bilin12_odd_E, ...
                               bilin13_even_E, bilin13_odd_E, ...
                               bilin23_even_E, bilin23_odd_E, ...
                               trilin_tau1_odd_E, trilin_tau2_even_E, ...
                               trilin_tau2_odd_E, trilin_tau3_even_E, ...
                               trilin_tau3_odd_E, trilin_tau4_even_E};

                    beads_C = {id_bead, bead1, bead2, bead3,  ...
                               bilin12_even_C, bilin12_odd_C, ...
                               bilin13_even_C, bilin13_odd_C, ...
                               bilin23_even_C, bilin23_odd_C, ...
                               trilin_tau1_odd_C, trilin_tau2_even_C, ...
                               trilin_tau2_odd_C, trilin_tau3_even_C, ...
                               trilin_tau3_odd_C, trilin_tau4_even_C};

                % Connected (entanglement-based) correlation function
                case "E"
                    beads = {id_bead, bead1, bead2, bead3, ...
                             bilin12_even_E, bilin12_odd_E, ...
                             bilin13_even_E, bilin13_odd_E, ...
                             bilin23_even_E, bilin23_odd_E, ...
                             trilin_tau1_odd_E, trilin_tau2_even_E, ...
                             trilin_tau2_odd_E, trilin_tau3_even_E, ...
                             trilin_tau3_odd_E, trilin_tau4_even_E};

                    beads_E = {blank, blank, blank, blank, ...
                               bilin12_even_E, bilin12_odd_E, ...
                               bilin13_even_E, bilin13_odd_E, ...
                               bilin23_even_E, bilin23_odd_E, ...
                               trilin_tau1_odd_E, trilin_tau2_even_E, ...
                               trilin_tau2_odd_E, trilin_tau3_even_E, ...
                               trilin_tau3_odd_E, trilin_tau4_even_E};

                    beads_C = {id_bead, bead1, bead2, bead3,  ...
                               blank, blank, ...
                               blank, blank, ...
                               blank, blank, ...
                               blank, blank, ...
                               blank, blank, ...
                               blank, blank};

                % Compound correlation function
                case "C"
                    beads = {id_bead, bead1, bead2, bead3, ...
                             bilin12_even_C, bilin12_odd_C, ...
                             bilin13_even_C, bilin13_odd_C, ...
                             bilin23_even_C, bilin23_odd_C, ...
                             trilin_tau1_odd_C, trilin_tau2_even_C, ...
                             trilin_tau2_odd_C, trilin_tau3_even_C, ...
                             trilin_tau3_odd_C, trilin_tau4_even_C};

                    beads_E = {blank, blank, blank, blank, ...
                               blank, blank, ...
                               blank, blank, ...
                               blank, blank, ...
                               blank, blank, ...
                               blank, blank, ...
                               blank, blank};

                    beads_C = {id_bead, bead1, bead2, bead3,  ...
                               bilin12_even_C, bilin12_odd_C, ...
                               bilin13_even_C, bilin13_odd_C, ...
                               bilin23_even_C, bilin23_odd_C, ...
                               trilin_tau1_odd_C, trilin_tau2_even_C, ...
                               trilin_tau2_odd_C, trilin_tau3_even_C, ...
                               trilin_tau3_odd_C, trilin_tau4_even_C};
            end

            % Set the axis handle of the i-th state to subplot i
            axes(ha(i));

            % Plot the individual beads
            for bb=1:size(beads,2)
                % Don't draw the id bead if not specified otherwise
                % Don't draw bilinear beads if they fall below the tolerance
                if (~Id && bb==1) || (bb > 3 && all(abs(beads{bb}(:)) < tolerance))
                    continue;
                end

                % Draw the bead surface plot
                s = surf(X+beads_positions(bb,1), Y+beads_positions(bb,2), ...
                    Z+beads_positions(bb,3), beads{bb});

                % Determine the color data and assign it to the surface plot
                set(s, 'CData', get_blended_colors(beads{bb}, beads_E{bb}, beads_C{bb}, GridPoints, false), ...
                    'EdgeColor','none');

                % Prevent perspective distortion
                axis equal;

                % The axis limits
                set(gca,'xlim', [-10, 8]);
                set(gca, 'ylim', [-7, 15]);
                set(gca, 'zlim', [-2, 2]);

                if viewpoint == "oblique" && EntArcs
                    set(gca,'xlim', [-11.5, 8]);
                end

                % Hide coordinate axes
                axis off;

                % The beads label
                if viewpoint ~= "none"
                    text(gca, label_positions(bb,1), label_positions(bb,2), label_positions(bb,3), ...
                        labels{bb}, 'Interpreter', 'Latex', 'FontSize', label_fontsize, ...
                        'VerticalAlignment','bottom');
                end

                hold on;
            end

            % Material properties
            shading interp;
            lightangle(100,45);
            material([1 0 0.2]);
            view(CameraAngle);

            hold on;

            % Draw the connector bars with radius proportional to the
            % entanglement of formation
            if EntArcs
                % Connector Q1-Q2
                % Y-symmetric Bézier curves (P0 and P3, and P1 and P2 have same
                % y-coordinates, respectively) have a maximum y_coordinate at
                % 3/4 y(P1) = 3/4 y(P2)

                if two_tangles(i,1) > tolerance
                    bezier(ha(i), {beads_positions(2,:), beads_positions(2,:) + [0, 4/3*beads_positions(5,2), 0], ...
                        beads_positions(3,:) + [0, 4/3*beads_positions(5,2), 0], beads_positions(3,:)}, ...
                        two_tangles(i, 1), [0 0 0]);
                end

                if two_tangles(i,2) > tolerance
                    % Connector Q1-Q3
                    bezier(ha(i), {beads_positions(2,:), beads_positions(2,:) + [0, 4/3*beads_positions(7,2), 0], ...
                        beads_positions(4,:) + [0, 4/3*beads_positions(7,2), 0], beads_positions(4,:)}, ...
                        two_tangles(i, 2), [0 0 0]);
                end

                if two_tangles(i,3) > tolerance
                    % Connector Q2-Q3
                    bezier(ha(i), {beads_positions(3,:), beads_positions(3,:) + [0, 4/3*beads_positions(9,2), 0], ...
                        beads_positions(4,:) + [0, 4/3*beads_positions(9,2), 0], beads_positions(4,:)}, ...
                        two_tangles(i, 3), [0 0 0]);
                end

                if (three_tangles(i) > tolerance) || (ent_norm_123(i) > tolerance)
                    tangle_color = [0 0 0];
                    ent_norm_abs = abs(three_tangles(i) - ent_norm_123(i));
                    if ent_norm_123(i) - three_tangles(i) > 0
                        norm_color = [.6 .35 .7];
                    else
                        norm_color = [.55, .6, .05];
                    end
                    if ent_norm_abs ~= 0
                        % Triple connectors Q1/Q2/Q3-τ1
                        % Set P2 to the same coordinates as P3 for Q1/3-τ1 connectors
                        bezier(ha(i), {beads_positions(2,:), beads_positions(2,:) + [0, beads_positions(11,2), 0], ...
                            beads_positions(11,:), beads_positions(11,:)}, ...
                            ent_norm_abs, norm_color, true, 5, true);

                        % Just a straight line for Q2-τ1
                        bezier(ha(i), {beads_positions(3,:), beads_positions(3,:), ...
                            beads_positions(11,:), beads_positions(11,:)}, ...
                            ent_norm_abs, norm_color, true, 5, true);

                        % Q3-τ1
                        bezier(ha(i), {beads_positions(4,:), beads_positions(4,:) + [0, beads_positions(11,2), 0], ...
                            beads_positions(11,:), beads_positions(11,:)}, ...
                            ent_norm_abs, norm_color, true, 5, true);
                    end

                    if three_tangles(i) ~= 0
                        % Triple connectors Q1/Q2/Q3-τ1
                        % Set P2 to the same coordinates as P3 for Q1/3-τ1 connectors
                        bezier(ha(i), {beads_positions(2,:), beads_positions(2,:) + [0, beads_positions(11,2), 0], ...
                            beads_positions(11,:), beads_positions(11,:)}, ...
                            three_tangles(i), tangle_color, true, 5, false);

                        % Just a straight line for Q2-τ1
                        bezier(ha(i), {beads_positions(3,:), beads_positions(3,:), ...
                            beads_positions(11,:), beads_positions(11,:)}, ...
                            three_tangles(i), tangle_color, true, 5, false);

                        % Q3-τ1
                        bezier(ha(i), {beads_positions(4,:), beads_positions(4,:) + [0, beads_positions(11,2), 0], ...
                            beads_positions(11,:), beads_positions(11,:)}, ...
                            three_tangles(i), tangle_color, true, 5, false);
                    end
                end
            end
            hold on;

            % Just draw the coordinate axis in the first subplot
            if i > 1
                continue;
            end

            % Calculate the cylinder and cone coordinates for every axis
            nTiles = 50;
            [x_cy, y_cy, z_cy] = cylinder(cylinder_radius, nTiles);
            cyl_pos = {x_cy, y_cy, z_cy};
            [x_cone, y_cone, z_cone] = cylinder([cone_radius, 0], nTiles);
            cone_pos = {x_cone, y_cone, z_cone*cone_height+1};
            disk_pos = {0, 0, 1};
            % The basic label position
            label_pos = {0,0,1.5};

            % Set gca to first subplot again
            axes(ha(1));
            for k = 1:3
                % Get the coordinates of this axis
                cyl_pos_temp = PermuteCoordinates(cyl_pos, axis_perm{k});
                cone_pos_temp = PermuteCoordinates(cone_pos, axis_perm{k});
                label_pos_temp = PermuteCoordinates(label_pos, axis_perm{k});
                disk_pos_temp = PermuteCoordinates(disk_pos, axis_perm{k});

                % The cylinder surface
                s = surf(gca, cyl_pos_temp{1}+axis_position(1), ...
                    cyl_pos_temp{2}+axis_position(2), ...
                    cyl_pos_temp{3}+axis_position(3));

                % Determine axis colors
                axcolor = zeros(size(x_cy,1), size(x_cy,2), 3);
                for rr = 1 : size(x_cy,1)
                    for cc = 1 : size(x_cy,2)
                        axcolor(rr, cc, :) = axis_colors(k, :);
                    end
                end

                % Assign colors to cylinder faces
                set(s, 'CData', axcolor,  'EdgeColor', 'none');

                hold on;

                % The cone surface
                s = surf(gca, cone_pos_temp{1}+axis_position(1), ...
                    cone_pos_temp{2}+axis_position(2), ...
                    cone_pos_temp{3}+axis_position(3));

                % Assign color
                set(s, 'CData', axcolor,  'EdgeColor', 'none');

                % Draw a disk at the base of the cone
                % The normal vector
                v = [label_pos_temp{1} label_pos_temp{2} label_pos_temp{3}];
                % Get the disk coordinates
                [xc, yc, zc] = get_filled_circle(v, cone_radius, nTiles);
                % Draw the disk
                surf(gca, disk_pos_temp{1}+axis_position(1)+xc, ...
                    disk_pos_temp{2}+axis_position(2)+yc, ...
                    disk_pos_temp{3}+axis_position(3)+zc, ...
                    'FaceColor', axis_colors(k, :), 'EdgeColor', 'none');

                % The axis label
                if viewpoint == "top" && k == 1
                    label_pos_temp = {-.4, 0, 1.5};
                elseif viewpoint == "top" && k == 2
                    label_pos_temp = {0, 1.7, 0};
                end
                text(label_pos_temp{1}+axis_position(1), ...
                    label_pos_temp{2}+axis_position(2), ...
                    label_pos_temp{3}+axis_position(3), axis_labels{k}, ...
                    'Color', axcolor(1,1,:), 'FontSize', label_fontsize, ...
                    'Parent', gca,'Interpreter','Latex', ...
                    'HorizontalAlignment','center');

                hold on;
            end
        end

        hold off

        % Hide last subplot axis if density operator has odd rank
        if mod(maxidx, 2)
            axes(ha(numel(ha)))
            axis off;
        end

        if maxidx < 8 && states > 8
            for aa = maxidx + 1 : 1 : 8
                axes(ha(aa))
                axis off;
            end
        end

        % Optionally draw lines between subplots
        if Lines
            % Draw vertical lines between columns
            for c = 1:nCols-1
                % Get the position of a subplot in this column (first row)
                pos = get(ha(c), 'Position');
                % X position of the right edge of this subplot
                x = pos(1) + pos(3) + gap(2)/2;
                if maxidx < 8 && states > 8
                    annotation(fig, 'line', [x x], [maxidx/8 1], 'Color', 'k', 'LineWidth', 0.5);
                else
                    annotation(fig, 'line', [x x], [0 1], 'Color', 'k', 'LineWidth', 0.5);
                end
            end

            % Draw horizontal lines between rows
            for r = 1:ceil(maxidx/nCols)-1
                % First column of this row
                idx = (r-1)*nCols + 1;
                pos = get(ha(idx), 'Position');
                y = pos(2) - gap(1)/2;
                annotation(fig, 'line', [0 1], [y y], 'Color', 'k', 'LineWidth', 0.5);
            end
        end

        % If desired, add text annotations in the top-left corner of each subplot
        if Weights
            for pp = 1:maxidx
                % Position of subplot handle
                sp_pos = get(ha(pp), 'Position');

                % Top-left corner (adjust x/y slightly to stay inside)
                x_shift = 0.01;
                y_shift = 0.055;
                if viewpoint == "top"
                    x_shift = 0;
                    y_shift = 0.055;
                end
                x = sp_pos(1) + x_shift;

                if nRows == 1
                    y = 1 - y_shift;
                else
                    y = round((nRows*nCols-pp+1)/2)*(1/nRows) - y_shift;
                end

                if Tangles
                    if viewpoint == "oblique"
                        value_string = {sprintf('$p_%d$ = %.2f\\%%, $\\tau_{12}$ = %.2f, $\\tau_{13}$ = %.2f, $\\tau_{23}$ = %.2f,',...
                            pp, probabilities(pp)*100, two_tangles(pp,1),two_tangles(pp,2),two_tangles(pp,3)),...
                            sprintf('$\\tau_{123}$ = %.2f',three_tangles(pp))};

                    else
                        value_string = sprintf(['$\\begin{array}{lcl}' ...
                            'p_%d & \\!\\!\\! = \\!\\!\\! & %.2f\\%% \\\\' ...
                            '\\tau_{12} & \\!\\!\\! = \\!\\!\\! & %.2f \\\\' ...
                            '\\tau_{13} & \\!\\!\\! = \\!\\!\\! & %.2f \\\\' ...
                            '\\tau_{23} & \\!\\!\\! = \\!\\!\\! & %.2f \\\\' ...
                            '\\tau_{123} & \\!\\!\\! = \\!\\!\\! & %.2f' ...
                            '\\end{array}$'], ...
                            pp, probabilities(pp)*100, ...
                            two_tangles(pp,1), two_tangles(pp,2), two_tangles(pp,3), ...
                            three_tangles(pp));
                    end

                else
                    value_string = sprintf('$p_%d$ = %.2f\\%%', pp, probabilities(pp)*100);
                end

                % Add the probability in the top-left corner
                annotation('textbox', [x y 0.05 0.05], ...
                    'String', value_string,...
                    'FitBoxToText', 'on', ...
                    'LineStyle', 'none', ...
                    'VerticalAlignment', 'top', ...
                    'HorizontalAlignment', 'left', ...
                    'FontSize', label_fontsize, ...
                    'Interpreter', 'latex');
            end
        end
        % Export graphic
        exportgraphics(gcf,path_,'Resolution',600);

        % Finally close the figure
        close;
        count = count + 8;
    end
end
