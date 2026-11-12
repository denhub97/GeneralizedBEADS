function [] = visualize3_universal(Beads, CameraAngle, BeadType, Id, ...
                                EntanglementArcs, TwoTangles, ThreeTangle, ...
                                EntNorm, GridPoints, RadialPlot, Normalize, Path, FontSize, ClassCorrRecursive)

    %!!! visualize3_universal visualizes the BEADS representation of a
    %    3-qubit state for one bead type with a spherical/radial plot. The
    %    absolute value of the correlation function corresponds to the color
    %    or the radius of the beads from the origin. This function only
    %    handles 3 qubit systems
    %
    %    inputs: Beads: the beads data
    %            CameraAngle: position of the camera / viewpoint
    %            BeadType: the type of correlation function beads (possible:
    %                      "T", "E", "C", "Sigma", "Sigma_Cl", "Sigma_D",
    %                      "D", "M", "M_Cl")
    %            Id: specify whether id bead shall be drawn (boolean)
    %            EntanglementArcs: specify whether connector arcs shall be
    %                              drawn (boolean)
    %            TwoTangles: the array of two-tangles
    %            ThreeTangle: the state's three tangle
    %            EntNorm: the state's tripartite entanglement norm
    %                     (to detect W-type entanglement)
    %            GridPoints: the beads resolution [theta, phi]
    %            RadialPlot: specify whether the beads shall be plotted with
    %                        absolute spherical function value as radius from
    %                        the origin
    %            Normalize: Normalize the droplets, this will not draw the
    %                       id droplet (boolean)
    %            Path: path to where the final image shall be saved
    %            FontSize: font size for all labels and annotations
    %            ClassCorrRecursive: (optional, default true) have the
    %                                projectors of the classical beads been
    %                                determined recursively for the marginals
    %                                and the tripartite structure (boolean),
    %                                this determines the target qubit labels
    %
    %    outputs: none, the figure is saved to Path

    % Handle default parameters
    if nargin < 14
        ClassCorrRecursive = true;
    end

    % Resolution of plot
    grid_points_theta = GridPoints(1);
    grid_points_phi = GridPoints(2);

    % Layout data
    label_shift = 0.7;
    label_fontsize = FontSize;
    tolerance = 1e-3;

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

    label_r = 1.3*sqrt(2*(1.2*0.7)^2);
    % Top view
    if CameraAngle(1) == 0 && CameraAngle(2) == 90
        viewpoint = "top";

        % Label positions
        label_positions = [label_shift*label_r,    0,                   0; % Id
                           label_shift*label_r,    0,                   0; % Q1
                           label_shift*label_r,    0,                   0;  % Q2
                           label_shift*label_r,    0,                   0;  % Q3
                           label_shift*cos(1.7*pi/4)*1.25, label_shift*sin(1.7*pi/4)*1.25,     0; % {1,2}_even
                           label_shift*cos(1.7*pi/4)*1.25, label_shift*sin(1.7*pi/4)*1.25,     0; % {1,2}_odd
                           label_shift*1.2,        label_shift*1.2,     0; % {1,3}_even
                           label_shift*1.2,        label_shift*1.2,     0; % {1,3}_odd
                           label_shift*cos(7.3*pi/8)*4.0-.2,1.2*label_shift*sin(7.3*pi/8)*4.0,0; % {2,3}_even
                           label_shift*cos(7.3*pi/8)*4.35,label_shift*sin(7.3*pi/8)*4.35,0; % {2,3}_odd
                           label_shift*1.2,        label_shift*1.2,     0; % Tau1_odd
                           label_shift*1.2,        label_shift*1.2,     0; % Tau2_even
                           label_shift*1.2,        label_shift*1.2,     0; % Tau2_odd
                           label_shift*1.2,        label_shift*1.2,     0; % Tau3_even
                           label_shift*1.2,        label_shift*1.2,     0; % Tau3_odd
                           label_shift*1.2,        label_shift*1.2,     0];% Tau4_even

    % Oblique views
    elseif CameraAngle(1) == 54 && CameraAngle(2) == 45
        viewpoint = "oblique";

        if EntanglementArcs
            % Label positions for the case that connector lines are enabled
            label_positions = [label_r*cos(pi/12),          label_r*sin(pi/12),        0; % Id
                               label_r*1.35*cos(pi/12),     label_r*1.35*sin(-pi/12),  0; % Q1
                               label_r*1.35*cos(pi/12),     label_r*1.35*sin(-pi/12),  0; % Q2
                               label_r*1.35*cos(pi/12),     label_r*1.35*sin(-pi/12),  0; % Q3
                               2.65*label_r*cos(11.7*pi/8), 2.65*label_r*sin(11.7*pi/8),0; % {1,2}_even
                               label_r*cos(7*pi/8),         label_r*sin(7*pi/8),       0; % {1,2}_odd
                               2.85*label_r*(-cos(pi/4)),   2.85*label_r*sin(-pi/4),   0; % {1,3}_even
                               1.5*label_r*cos(pi/12),      1.5*label_r*sin(-pi/12),   0; % {1,3}_odd
                               2.65*label_r*cos(11.7*pi/8), 2.65*label_r*sin(11.7*pi/8),0; % {2,3}_even
                               1.35*label_r*(cos(4.25*pi/4)),1.35*label_r*sin(4.25*pi/4),     0; % {2,3}_odd
                               label_r*cos(pi/12),          label_r*sin(pi/12),        0; % Tau1_odd
                               label_r*cos(pi/12),          label_r*sin(pi/12),        0; % Tau2_even
                               label_r*cos(pi/12),          label_r*sin(pi/12),        0; % Tau2_odd
                               label_r*cos(pi/12),          label_r*sin(pi/12),        0; % Tau3_even
                               label_r*cos(pi/12),          label_r*sin(pi/12),        0; % Tau3_odd
                               label_r*cos(pi/12),          label_r*sin(pi/12),        0];% Tau4_even
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

    % Assign data to the individual beads
    id_bead = Beads{1};
    bead1 = Beads{2};
    bead2 = Beads{3};
    bead3 = Beads{4};

    % Correlation function types for even and odd case in the following order:
    % T-Bead: Total correlation function
    % E-Bead: Connected correlation function (quantum entanglement)
    % C-Bead: Compound, C = T - E
    % Sigma-Bead: Compound - reducible components of n-partite correlation
    %          functions (e.g. components which factorize into components of lower
    %          linearity)
    % SigmaCl-Bead: Separable contribution to classical correlation function.
    %          Projectors are optimized such that measurement on target qubit
    %          minimally disturbs the overall system
    %          (in other words: such that information on remaining subsystem is maximal)
    % SigmaD-Bead: Non-entanglement-based (quantum) correlation function
    %         which is nonclassical
    % M-Bead: (Total) mutual information: total correlation function where reducible
    %         components of n-partite correlation functions have been subtracted
    % MCl-Bead: Classical mutual information correlation function. M-Beads are projected
    %          using optimized projectors (same as for SigmaCl-Beads)
    % ECl-Bead: Connected correlation function projected using the optimized
    %          projectors (same as for SigmaCl-Beads)

    bilin12_even_T = Beads{5};
    bilin12_even_E = Beads{6};
    bilin12_even_C = Beads{7};
    bilin12_even_Sigma = Beads{8};
    bilin12_even_SigmaCl = Beads{9};
    bilin12_even_SigmaD = Beads{10};
    bilin12_even_M = Beads{11};
    bilin12_even_MCl = Beads{12};

    bilin13_even_T = Beads{13};
    bilin13_even_E = Beads{14};
    bilin13_even_C = Beads{15};
    bilin13_even_Sigma = Beads{16};
    bilin13_even_SigmaCl = Beads{17};
    bilin13_even_SigmaD = Beads{18};
    bilin13_even_M = Beads{19};
    bilin13_even_MCl = Beads{20};

    bilin23_even_T = Beads{21};
    bilin23_even_E = Beads{22};
    bilin23_even_C = Beads{23};
    bilin23_even_Sigma = Beads{24};
    bilin23_even_SigmaCl = Beads{25};
    bilin23_even_SigmaD = Beads{26};
    bilin23_even_M = Beads{27};
    bilin23_even_MCl = Beads{28};

    trilin_tau2_even_T = Beads{29};
    trilin_tau2_even_E = Beads{30};
    trilin_tau2_even_C = Beads{31};
    trilin_tau2_even_Sigma = Beads{32};
    trilin_tau2_even_SigmaCl = Beads{33};
    trilin_tau2_even_SigmaD = Beads{34};
    trilin_tau2_even_M = Beads{35};
    trilin_tau2_even_MCl = Beads{36};

    trilin_tau3_even_T = Beads{37};
    trilin_tau3_even_E = Beads{38};
    trilin_tau3_even_C = Beads{39};
    trilin_tau3_even_Sigma = Beads{40};
    trilin_tau3_even_SigmaCl = Beads{41};
    trilin_tau3_even_SigmaD = Beads{42};
    trilin_tau3_even_M = Beads{43};
    trilin_tau3_even_MCl = Beads{44};

    trilin_tau4_even_T = Beads{45};
    trilin_tau4_even_E = Beads{46};
    trilin_tau4_even_C = Beads{47};
    trilin_tau4_even_Sigma = Beads{48};
    trilin_tau4_even_SigmaCl = Beads{49};
    trilin_tau4_even_SigmaD = Beads{50};
    trilin_tau4_even_M = Beads{51};
    trilin_tau4_even_MCl = Beads{52};

    bilin12_odd_T = Beads{53};
    bilin12_odd_E = Beads{54};
    bilin12_odd_C = Beads{55};
    bilin12_odd_Sigma = Beads{56};
    bilin12_odd_SigmaCl = Beads{57};
    bilin12_odd_SigmaD = Beads{58};
    bilin12_odd_M = Beads{59};
    bilin12_odd_MCl = Beads{60};

    bilin13_odd_T = Beads{61};
    bilin13_odd_E = Beads{62};
    bilin13_odd_C = Beads{63};
    bilin13_odd_Sigma = Beads{64};
    bilin13_odd_SigmaCl = Beads{65};
    bilin13_odd_SigmaD = Beads{66};
    bilin13_odd_M = Beads{67};
    bilin13_odd_MCl = Beads{68};

    bilin23_odd_T = Beads{69};
    bilin23_odd_E = Beads{70};
    bilin23_odd_C = Beads{71};
    bilin23_odd_Sigma = Beads{72};
    bilin23_odd_SigmaCl = Beads{73};
    bilin23_odd_SigmaD = Beads{74};
    bilin23_odd_M = Beads{75};
    bilin23_odd_MCl = Beads{76};

    trilin_tau2_odd_T = Beads{77};
    trilin_tau2_odd_E = Beads{78};
    trilin_tau2_odd_C = Beads{79};
    trilin_tau2_odd_Sigma = Beads{80};
    trilin_tau2_odd_SigmaCl = Beads{81};
    trilin_tau2_odd_SigmaD = Beads{82};
    trilin_tau2_odd_M = Beads{83};
    trilin_tau2_odd_MCl = Beads{84};

    trilin_tau3_odd_T = Beads{85};
    trilin_tau3_odd_E = Beads{86};
    trilin_tau3_odd_C = Beads{87};
    trilin_tau3_odd_Sigma = Beads{88};
    trilin_tau3_odd_SigmaCl = Beads{89};
    trilin_tau3_odd_SigmaD = Beads{90};
    trilin_tau3_odd_M = Beads{91};
    trilin_tau3_odd_MCl = Beads{92};

    trilin_tau1_odd_T = Beads{93};
    trilin_tau1_odd_E = Beads{94};
    trilin_tau1_odd_C = Beads{95};
    trilin_tau1_odd_Sigma = Beads{96};
    trilin_tau1_odd_SigmaCl = Beads{97};
    trilin_tau1_odd_SigmaD = Beads{98};
    trilin_tau1_odd_M = Beads{99};
    trilin_tau1_odd_MCl = Beads{100};

    bilin12_even_ECl = Beads{101};
    bilin12_odd_ECl = Beads{102};
    bilin13_even_ECl = Beads{103};
    bilin13_odd_ECl = Beads{104};
    bilin23_even_ECl = Beads{105};
    bilin23_odd_ECl = Beads{106};
    trilin_tau1_odd_ECl = Beads{107};
    trilin_tau2_odd_ECl = Beads{108};
    trilin_tau2_even_ECl = Beads{109};
    trilin_tau3_odd_ECl = Beads{110};
    trilin_tau3_even_ECl = Beads{111};
    trilin_tau4_even_ECl = Beads{112};

    % The target qubits of the classical beads, [q12, q13, q23, q1, q2]
    targetqubits = Beads{113};

    % Labels of the target qubits for the beads labels
    targetlabel12 = num2str(targetqubits(1));
    targetlabel13 = num2str(targetqubits(2));
    targetlabel23 = num2str(targetqubits(3));
    targetlabel123 = append(num2str(targetqubits(4)),',',num2str(targetqubits(5)));

    % If the projectors have not been determined recursively, the same target
    % qubits apply to all beads
    if ~ClassCorrRecursive
        tq = targetqubits(4:5);
        if tq(1) == tq(2)
            tq = tq(1);
            targetlabel12 = num2str(tq);
            targetlabel13 = num2str(tq);
            targetlabel23 = num2str(tq);
            targetlabel123 = append(num2str(tq));
        else
            targetlabel12 = append(num2str(tq(1)),',',num2str(tq(2)));
            targetlabel13 = append(num2str(tq(1)),',',num2str(tq(2)));
            targetlabel23 = append(num2str(tq(1)),',',num2str(tq(2)));
            targetlabel123 = append(num2str(tq(1)),',',num2str(tq(2)));
        end
    end

    % Determine the beads labels
    switch BeadType
        case "T"
            labels = {'Id','$Q_1$','$Q_2$','$Q_3$','$T\{1,2\}_{even}$',...
                      '$T\{1,2\}_{odd}$','$T\{1,3\}_{even}$','$T\{1,3\}_{odd}$',...
                      '$T\{2,3\}_{even}$','$T\{2,3\}_{odd}$',...
                      '$T\{\tau_1\}_{odd}$','$T\{\tau_2\}_{even}$',...
                      '$T\{\tau_2\}_{odd}$','$T\{\tau_3\}_{even}$',...
                      '$T\{\tau_3\}_{odd}$','$T\{\tau_4\}_{even}$'};

        case "E"
            labels = {'Id','$Q_1$','$Q_2$','$Q_3$','$E\{1,2\}_{even}$',...
                      '$E\{1,2\}_{odd}$','$E\{1,3\}_{even}$','$E\{1,3\}_{odd}$',...
                      '$E\{2,3\}_{even}$','$E\{2,3\}_{odd}$',...
                      '$E\{\tau_1\}_{odd}$','$E\{\tau_2\}_{even}$',...
                      '$E\{\tau_2\}_{odd}$','$E\{\tau_3\}_{even}$',...
                      '$E\{\tau_3\}_{odd}$','$E\{\tau_4\}_{even}$'};

        case "C"
            labels = {'Id','$Q_1$','$Q_2$','$Q_3$','$C\{1,2\}_{even}$',...
                      '$C\{1,2\}_{odd}$','$C\{1,3\}_{even}$','$C\{1,3\}_{odd}$',...
                      '$C\{2,3\}_{even}$','$C\{2,3\}_{odd}$',...
                      '$C\{\tau_1\}_{odd}$','$C\{\tau_2\}_{even}$',...
                      '$C\{\tau_2\}_{odd}$','$C\{\tau_3\}_{even}$',...
                      '$C\{\tau_3\}_{odd}$','$C\{\tau_4\}_{even}$'};

        case "Sigma"
            labels = {'Id','$Q_1$','$Q_2$','$Q_3$','$\Sigma\{1,2\}_{even}$',...
                      '$\Sigma\{1,2\}_{odd}$','$\Sigma\{1,3\}_{even}$','$\Sigma\{1,3\}_{odd}$',...
                      '$\Sigma\{2,3\}_{even}$','$\Sigma\{2,3\}_{odd}$',...
                      '$\Sigma\{\tau_1\}_{odd}$','$\Sigma\{\tau_2\}_{even}$',...
                      '$\Sigma\{\tau_2\}_{odd}$','$\Sigma\{\tau_3\}_{even}$',...
                      '$\Sigma\{\tau_3\}_{odd}$','$\Sigma\{\tau_4\}_{even}$'};

        case "Sigma_Cl"
            labels = {'Id','$Q_1$','$Q_2$','$Q_3$', ...
                      append('$\Sigma_{Cl}^{(',targetlabel12,')}\{1,2\}_{even}$'), ...
                      append('$\Sigma_{Cl}^{(',targetlabel12,')}\{1,2\}_{odd}$'), ...
                      append('$\Sigma_{Cl}^{(',targetlabel13,')}\{1,3\}_{even}$'), ...
                      append('$\Sigma_{Cl}^{(',targetlabel13,')}\{1,3\}_{odd}$'),...
                      append('$\Sigma_{Cl}^{(',targetlabel23,')}\{2,3\}_{even}$'), ...
                      append('$\Sigma_{Cl}^{(',targetlabel23,')}\{2,3\}_{odd}$'),...
                      append('$\Sigma_{Cl}^{(',targetlabel123,')}\{1,2,3\ \tau_1\}_{odd}$'), ...
                      append('$\Sigma_{Cl}^{(',targetlabel123,')}\{1,2,3\ \tau_2\}_{even}$'), ...
                      append('$\Sigma_{Cl}^{(',targetlabel123,')}\{1,2,3\ \tau_2\}_{odd}$'), ...
                      append('$\Sigma_{Cl}^{(',targetlabel123,')}\{1,2,3\ \tau_3\}_{even}$'), ...
                      append('$\Sigma_{Cl}^{(',targetlabel123,')}\{1,2,3\ \tau_3\}_{odd}$'), ...
                      append('$\Sigma_{Cl}^{(',targetlabel123,')}\{1,2,3\ \tau_4\}_{even}$')};

        case "Sigma_D"
            labels = {'Id','$Q_1$','$Q_2$','$Q_3$', ...
                      append('$\Sigma_{D}^{(',targetlabel12,')}\{1,2\}_{even}$'), ...
                      append('$\Sigma_{D}^{(',targetlabel12,')}\{1,2\}_{odd}$'), ...
                      append('$\Sigma_{D}^{(',targetlabel13,')}\{1,3\}_{even}$'), ...
                      append('$\Sigma_{D}^{(',targetlabel13,')}\{1,3\}_{odd}$'),...
                      append('$\Sigma_{D}^{(',targetlabel23,')}\{2,3\}_{even}$'), ...
                      append('$\Sigma_{D}^{(',targetlabel23,')}\{2,3\}_{odd}$'),...
                      append('$\Sigma_{D}^{(',targetlabel123,')}\{1,2,3\ \tau_1\}_{odd}$'), ...
                      append('$\Sigma_{D}^{(',targetlabel123,')}\{1,2,3\ \tau_2\}_{even}$'), ...
                      append('$\Sigma_{D}^{(',targetlabel123,')}\{1,2,3\ \tau_2\}_{odd}$'), ...
                      append('$\Sigma_{D}^{(',targetlabel123,')}\{1,2,3\ \tau_3\}_{even}$'), ...
                      append('$\Sigma_{D}^{(',targetlabel123,')}\{1,2,3\ \tau_3\}_{odd}$'), ...
                      append('$\Sigma_{D}^{(',targetlabel123,')}\{1,2,3\ \tau_4\}_{even}$')};

        case "D"
            labels = {'Id','$Q_1$','$Q_2$','$Q_3$', ...
                      append('$D^{(',targetlabel12,')}\{1,2\}_{even}$'), ...
                      append('$D^{(',targetlabel12,')}\{1,2\}_{odd}$'), ...
                      append('$D^{(',targetlabel13,')}\{1,3\}_{even}$'), ...
                      append('$D^{(',targetlabel13,')}\{1,3\}_{odd}$'),...
                      append('$D^{(',targetlabel23,')}\{2,3\}_{even}$'), ...
                      append('$D^{(',targetlabel23,')}\{2,3\}_{odd}$'),...
                      append('$D^{(',targetlabel123,')}\{1,2,3\ \tau_1\}_{odd}$'), ...
                      append('$D^{(',targetlabel123,')}\{1,2,3\ \tau_2\}_{even}$'), ...
                      append('$D^{(',targetlabel123,')}\{1,2,3\ \tau_2\}_{odd}$'), ...
                      append('$D^{(',targetlabel123,')}\{1,2,3\ \tau_3\}_{even}$'), ...
                      append('$D^{(',targetlabel123,')}\{1,2,3\ \tau_3\}_{odd}$'), ...
                      append('$D^{(',targetlabel123,')}\{1,2,3\ \tau_4\}_{even}$')};

        case "M"
            labels = {'Id','$Q_1$','$Q_2$','$Q_3$','$M\{1,2\}_{even}$',...
                      '$M\{1,2\}_{odd}$','$M\{1,3\}_{even}$','$M\{1,3\}_{odd}$',...
                      '$M\{2,3\}_{even}$','$M\{2,3\}_{odd}$',...
                      '$M\{\tau_1\}_{odd}$','$M\{\tau_2\}_{even}$',...
                      '$M\{\tau_2\}_{odd}$','$M\{\tau_3\}_{even}$',...
                      '$M\{\tau_3\}_{odd}$','$M\{\tau_4\}_{even}$'};

        case "M_Cl"
            labels = {'Id','$Q_1$','$Q_2$','$Q_3$', ...
                      append('$M_{Cl}^{(',targetlabel12,')}\{1,2\}_{even}$'), ...
                      append('$M_{Cl}^{(',targetlabel12,')}\{1,2\}_{odd}$'), ...
                      append('$M_{Cl}^{(',targetlabel13,')}\{1,3\}_{even}$'), ...
                      append('$M_{Cl}^{(',targetlabel13,')}\{1,3\}_{odd}$'),...
                      append('$M_{Cl}^{(',targetlabel23,')}\{2,3\}_{even}$'), ...
                      append('$M_{Cl}^{(',targetlabel23,')}\{2,3\}_{odd}$'),...
                      append('$M_{Cl}^{(',targetlabel123,')}\{1,2,3\ \tau_1\}_{odd}$'), ...
                      append('$M_{Cl}^{(',targetlabel123,')}\{1,2,3\ \tau_2\}_{even}$'), ...
                      append('$M_{Cl}^{(',targetlabel123,')}\{1,2,3\ \tau_2\}_{odd}$'), ...
                      append('$M_{Cl}^{(',targetlabel123,')}\{1,2,3\ \tau_3\}_{even}$'), ...
                      append('$M_{Cl}^{(',targetlabel123,')}\{1,2,3\ \tau_3\}_{odd}$'), ...
                      append('$M_{Cl}^{(',targetlabel123,')}\{1,2,3\ \tau_4\}_{even}$')};

    end

    % Coordinate axes layout parameters
    cylinder_radius = 0.05;
    cone_height = 0.2;
    cone_radius = 0.1;
    axis_position = [-9.5, 0, 0];
    if viewpoint == "oblique" && EntanglementArcs
        axis_position = [-11, 0, 0];
    end

    % Permutations to get all 3D coordinate axes
    axis_perm = {[1 2 3], [1 3 2], [3 2 1]};
    axis_colors = [0  0  0;
                  .9 .8  0;
                   1  0  0];
    axis_labels = {'z', 'y', 'x'};

    % Define spherical coordinate mesh for plotting
    theta = linspace(0,pi,grid_points_theta);
    phi = linspace(0,2*pi,grid_points_phi);
    [Phi,Theta] = meshgrid(phi,theta);
    radius=1.0;

    % Cartesian from spherical coordinates for 3D meshgrid
    X = cos(Phi).*sin(Theta).*radius;
    Y = sin(Phi).*sin(Theta).*radius;
    Z = cos(Theta)*radius;

    % All beads in the same plot
    fig = figure;
    fig.Position = [0,0,800,600];
    % Set background color to white
    fig.Color = 'w';

    % One axis object
    ax = axes();
    axes(ax);
    % Facilitate plotting into ax with external functions (bezier)
    set(ax, 'NextPlot', 'add');

    % Zero array, required for color assignment
    blank = 0*X;

    % Determine the beads' spherical function values (beads) and values
    % used for determination of color (beads_E/C)
    switch BeadType
        % Total correlation function
        case "T"
            beads = {id_bead, bead1, bead2, bead3, ...
                     bilin12_even_T, bilin12_odd_T, ...
                     bilin13_even_T, bilin13_odd_T, ...
                     bilin23_even_T, bilin23_odd_T, ...
                     trilin_tau1_odd_T, trilin_tau2_even_T, ...
                     trilin_tau2_odd_T, trilin_tau3_even_T, ...
                     trilin_tau3_odd_T, trilin_tau4_even_T};

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

        % Correlation functions of separable components of the state
        % (classical and nonclassical correlations)
        case "Sigma"
            beads = {id_bead, bead1, bead2, bead3, ...
                     bilin12_even_Sigma, bilin12_odd_Sigma, ...
                     bilin13_even_Sigma, bilin13_odd_Sigma, ...
                     bilin23_even_Sigma, bilin23_odd_Sigma, ...
                     trilin_tau1_odd_Sigma, trilin_tau2_even_Sigma, ...
                     trilin_tau2_odd_Sigma, trilin_tau3_even_Sigma, ...
                     trilin_tau3_odd_Sigma, trilin_tau4_even_Sigma};

            beads_E = {blank, blank, blank, blank, ...
                       blank, blank, ...
                       blank, blank, ...
                       blank, blank, ...
                       blank, blank, ...
                       blank, blank, ...
                       blank, blank};

            beads_C = {id_bead, bead1, bead2, bead3,  ...
                       bilin12_even_Sigma, bilin12_odd_Sigma, ...
                       bilin13_even_Sigma, bilin13_odd_Sigma, ...
                       bilin23_even_Sigma, bilin23_odd_Sigma, ...
                       trilin_tau1_odd_Sigma, trilin_tau2_even_Sigma, ...
                       trilin_tau2_odd_Sigma, trilin_tau3_even_Sigma, ...
                       trilin_tau3_odd_Sigma, trilin_tau4_even_Sigma};

        % Correlation functions of separable components of the state
        % contributing to classical correlation (CMI)
        case "Sigma_Cl"
            beads = {id_bead, bead1, bead2, bead3, ...
                     bilin12_even_SigmaCl, bilin12_odd_SigmaCl, ...
                     bilin13_even_SigmaCl, bilin13_odd_SigmaCl, ...
                     bilin23_even_SigmaCl, bilin23_odd_SigmaCl, ...
                     trilin_tau1_odd_SigmaCl, trilin_tau2_even_SigmaCl, ...
                     trilin_tau2_odd_SigmaCl, trilin_tau3_even_SigmaCl, ...
                     trilin_tau3_odd_SigmaCl, trilin_tau4_even_SigmaCl};

            beads_E = {blank, blank, blank, blank, ...
                       blank, blank, ...
                       blank, blank, ...
                       blank, blank, ...
                       blank, blank, ...
                       blank, blank, ...
                       blank, blank};

            beads_C = {id_bead, bead1, bead2, bead3,  ...
                       bilin12_even_SigmaCl, bilin12_odd_SigmaCl, ...
                       bilin13_even_SigmaCl, bilin13_odd_SigmaCl, ...
                       bilin23_even_SigmaCl, bilin23_odd_SigmaCl, ...
                       trilin_tau1_odd_SigmaCl, trilin_tau2_even_SigmaCl, ...
                       trilin_tau2_odd_SigmaCl, trilin_tau3_even_SigmaCl, ...
                       trilin_tau3_odd_SigmaCl, trilin_tau4_even_SigmaCl};

        % Nonclassical correlation functions of separable components of the
        % state
        case "D"
            beads = {id_bead, bead1, bead2, bead3, ...
                     bilin12_even_M-bilin12_even_MCl, bilin12_odd_M-bilin12_odd_MCl, ...
                     bilin13_even_M-bilin13_even_MCl, bilin13_odd_M-bilin13_odd_MCl, ...
                     bilin23_even_M-bilin23_even_MCl, bilin23_odd_M-bilin23_odd_MCl, ...
                     trilin_tau1_odd_M-trilin_tau1_odd_MCl, trilin_tau2_even_M-trilin_tau2_even_MCl, ...
                     trilin_tau2_odd_M-trilin_tau2_odd_MCl, trilin_tau3_even_M-trilin_tau3_even_MCl, ...
                     trilin_tau3_odd_M-trilin_tau3_odd_MCl, trilin_tau4_even_M-trilin_tau4_even_MCl};

            beads_E = {blank, blank, blank, blank, ...
                       bilin12_even_E-bilin12_even_ECl, bilin12_odd_E-bilin12_odd_ECl, ...
                       bilin13_even_E-bilin13_even_ECl,  bilin13_odd_E-bilin13_odd_ECl, ...
                       bilin23_even_E-bilin23_even_ECl, bilin23_odd_E-bilin23_odd_ECl, ...
                       trilin_tau1_odd_E-trilin_tau1_odd_ECl, trilin_tau2_even_E-trilin_tau2_even_ECl, ...
                       trilin_tau2_odd_E-trilin_tau2_odd_ECl, trilin_tau3_even_E-trilin_tau3_even_ECl, ...
                       trilin_tau3_odd_E-trilin_tau3_odd_ECl, trilin_tau4_even_E-trilin_tau4_even_ECl};

            beads_C = {id_bead, bead1, bead2, bead3,  ...
                       bilin12_even_SigmaD, bilin12_odd_SigmaD, ...
                       bilin13_even_SigmaD, bilin13_odd_SigmaD, ...
                       bilin23_even_SigmaD, bilin23_odd_SigmaD, ...
                       trilin_tau1_odd_SigmaD, trilin_tau2_even_SigmaD, ...
                       trilin_tau2_odd_SigmaD, trilin_tau3_even_SigmaD, ...
                       trilin_tau3_odd_SigmaD, trilin_tau4_even_SigmaD};

        case "Sigma_D"
            beads = {id_bead, bead1, bead2, bead3, ...
                     bilin12_even_SigmaD, bilin12_odd_SigmaD, ...
                     bilin13_even_SigmaD, bilin13_odd_SigmaD, ...
                     bilin23_even_SigmaD, bilin23_odd_SigmaD, ...
                     trilin_tau1_odd_SigmaD, trilin_tau2_even_SigmaD, ...
                     trilin_tau2_odd_SigmaD, trilin_tau3_even_SigmaD, ...
                     trilin_tau3_odd_SigmaD, trilin_tau4_even_SigmaD};

            beads_E = {blank, blank, blank, blank, ...
                       blank, blank, ...
                       blank, blank, ...
                       blank, blank, ...
                       blank, blank, ...
                       blank, blank, ...
                       blank, blank};

            beads_C = {id_bead, bead1, bead2, bead3,  ...
                       bilin12_even_SigmaD, bilin12_odd_SigmaD, ...
                       bilin13_even_SigmaD, bilin13_odd_SigmaD, ...
                       bilin23_even_SigmaD, bilin23_odd_SigmaD, ...
                       trilin_tau1_odd_SigmaD, trilin_tau2_even_SigmaD, ...
                       trilin_tau2_odd_SigmaD, trilin_tau3_even_SigmaD, ...
                       trilin_tau3_odd_SigmaD, trilin_tau4_even_SigmaD};

        % Classical mutual information function
        case "M_Cl"
            beads = {id_bead, bead1, bead2, bead3, ...
                     bilin12_even_MCl, bilin12_odd_MCl, ...
                     bilin13_even_MCl, bilin13_odd_MCl, ...
                     bilin23_even_MCl, bilin23_odd_MCl, ...
                     trilin_tau1_odd_MCl, trilin_tau2_even_MCl, ...
                     trilin_tau2_odd_MCl, trilin_tau3_even_MCl, ...
                     trilin_tau3_odd_MCl, trilin_tau4_even_MCl};

            beads_E = {blank, blank, blank, blank, ...
                       bilin12_even_ECl, bilin12_odd_ECl, ...
                       bilin13_even_ECl, bilin13_odd_ECl, ...
                       bilin23_even_ECl, bilin23_odd_ECl, ...
                       trilin_tau1_odd_ECl, trilin_tau2_even_ECl, ...
                       trilin_tau2_odd_ECl, trilin_tau3_even_ECl, ...
                       trilin_tau3_odd_ECl, trilin_tau4_even_ECl};

            beads_C = {id_bead, bead1, bead2, bead3,  ...
                       bilin12_even_SigmaCl, bilin12_odd_SigmaCl, ...
                       bilin13_even_SigmaCl, bilin13_odd_SigmaCl, ...
                       bilin23_even_SigmaCl, bilin23_odd_SigmaCl, ...
                       trilin_tau1_odd_SigmaCl, trilin_tau2_even_SigmaCl, ...
                       trilin_tau2_odd_SigmaCl, trilin_tau3_even_SigmaCl, ...
                       trilin_tau3_odd_SigmaCl, trilin_tau4_even_SigmaCl};

        % Total mutual information correlation function
        case "M"
            beads = {id_bead, bead1, bead2, bead3, ...
                     bilin12_even_M, bilin12_odd_M, ...
                     bilin13_even_M, bilin13_odd_M, ...
                     bilin23_even_M, bilin23_odd_M, ...
                     trilin_tau1_odd_M, trilin_tau2_even_M, ...
                     trilin_tau2_odd_M, trilin_tau3_even_M, ...
                     trilin_tau3_odd_M, trilin_tau4_even_M};

            beads_E = {blank, blank, blank, blank, ...
                       bilin12_even_E, bilin12_odd_E, ...
                       bilin13_even_E, bilin13_odd_E, ...
                       bilin23_even_E, bilin23_odd_E, ...
                       trilin_tau1_odd_E, trilin_tau2_even_E, ...
                       trilin_tau2_odd_E, trilin_tau3_even_E, ...
                       trilin_tau3_odd_E, trilin_tau4_even_E};

            beads_C = {id_bead, bead1, bead2, bead3,  ...
                       bilin12_even_Sigma, bilin12_odd_Sigma, ...
                       bilin13_even_Sigma, bilin13_odd_Sigma, ...
                       bilin23_even_Sigma, bilin23_odd_Sigma, ...
                       trilin_tau1_odd_Sigma, trilin_tau2_even_Sigma, ...
                       trilin_tau2_odd_Sigma, trilin_tau3_even_Sigma, ...
                       trilin_tau3_odd_Sigma, trilin_tau4_even_Sigma};

    end

    % Normalize the beads, i.e., scale the radius (radial plot) or the color
    % data (otherwise) by the automatically determined scale factor
    radius_scale = 1;
    if Normalize
        radius_scale = autoscale(beads);
        if radius_scale ~= 1.0
            annotation("textbox",[0.273, 0.71, 0.07, 0.0425], 'String', append('$\times$', sprintf('%.2f', radius_scale)), 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
        end
    end

    % Plot the individual beads
    for bb = 1:size(beads,2)
        % The beads label
        if viewpoint ~= "none" && ...
           ((bb > 4 && any(abs(beads{bb}(:)) > tolerance)) || ...
           (bb > 1 && bb <= 4) || (bb == 1 && Id))
            text_fontsize = label_fontsize;
            % Draw the bead label
            text(ax, beads_positions(bb,1)+label_positions(bb,1), ...
                     beads_positions(bb,2)+label_positions(bb,2), ...
                     beads_positions(bb,3)+label_positions(bb,3), ...
                     labels{bb}, 'Interpreter', 'Latex', 'FontSize', text_fontsize, ...
                     'VerticalAlignment','bottom', 'HorizontalAlignment', 'left');
        end

        % Don't draw the id bead if not specified otherwise
        % Don't draw bilinear beads if they fall below the tolerance
        if ((~Id || Normalize) && bb == 1) || ...
           (bb > 4 && all(abs(beads{bb}(:)) < tolerance))
            continue;
        end

        % The radii array
        if RadialPlot
            r = abs(beads{bb});
        else
            r = ones(grid_points_theta, grid_points_phi);
        end
        % Draw the bead surface plot
        if RadialPlot
            s = surf(ax, radius_scale*r.*X+beads_positions(bb,1), radius_scale*r.*Y+beads_positions(bb,2), ...
                         radius_scale*r.*Z+beads_positions(bb,3), beads{bb});
        else
            s = surf(ax, r.*X+beads_positions(bb,1), r.*Y+beads_positions(bb,2), ...
                         r.*Z+beads_positions(bb,3), radius_scale*beads{bb});
        end

        % Determine the color data and assign it to the surface plot
        set(s, 'CData', get_blended_colors(radius_scale*beads{bb}, radius_scale*beads_E{bb}, radius_scale*beads_C{bb}, GridPoints, RadialPlot), ...
            'EdgeColor','none');

        hold on;
    end

    % Material properties
    shading interp;
    lightangle(100,45);
    if RadialPlot
        material ([0.7 0.6 0.5]);
    else
        material([1 0 0.2]);
    end
    view(CameraAngle);
    axis off
    axis equal
    % The axis limits
    set(gca,'xlim', [-10, 8]);
    set(gca, 'ylim', [-7, 15]);
    set(gca, 'zlim', [-2, 2]);
    if viewpoint == "oblique" && EntanglementArcs
          set(gca,'xlim', [-11.5, 8]);
    end

    % Draw the connector bars with radius proportional to the
    % entanglement of formation
    if EntanglementArcs
        % Connector Q1-Q2
        % Y-symmetric Bézier curves (P0 and P3, and P1 and P2 have same
        % y-coordinates, respectively) have a maximum y_coordinate at
        % 3/4 y(P1) = 3/4 y(P2)
        if TwoTangles(1) > tolerance
            % Connector Q1-Q2
            bezier(ax, {beads_positions(2,:), beads_positions(2,:) + [0, 4/3*beads_positions(5,2), 0], ...
                        beads_positions(3,:) + [0, 4/3*beads_positions(5,2), 0], beads_positions(3,:)}, ...
                        TwoTangles(1), [0 0 0]);
        end

        if TwoTangles(2) > tolerance
            % Connector Q1-Q3
            bezier(ax, {beads_positions(2,:), beads_positions(2,:) + [0, 4/3*beads_positions(7,2), 0], ...
                        beads_positions(4,:) + [0, 4/3*beads_positions(7,2), 0], beads_positions(4,:)}, ...
                        TwoTangles(2), [0 0 0]);
        end

        if TwoTangles(3) > tolerance
            % Connector Q2-Q3
            bezier(ax, {beads_positions(3,:), beads_positions(3,:) + [0, 4/3*beads_positions(9,2), 0], ...
                        beads_positions(4,:) + [0, 4/3*beads_positions(9,2), 0], beads_positions(4,:)}, ...
                        TwoTangles(3), [0 0 0]);

        end

        if (ThreeTangle > tolerance) || (EntNorm > tolerance)
            tangle_color = [0 0 0];
            ent_norm_abs = abs(ThreeTangle - EntNorm);
            if EntNorm - ThreeTangle > 0
                norm_color = [.6 .35 .7];
            else
                norm_color = [.55, .6, .05];
            end
            if ent_norm_abs ~= 0
                % Triple connectors Q1/Q2/Q3-τ1
                % Set P2 to the same coordinates as P3 for Q1/3-τ1 connectors
                bezier(ax, {beads_positions(2,:), beads_positions(2,:) + [0, beads_positions(11,2), 0], ...
                            beads_positions(11,:), beads_positions(11,:)}, ...
                            ent_norm_abs, norm_color, true, 5, true);

                % Just a straight line for Q2-τ1
                bezier(ax, {beads_positions(3,:), beads_positions(3,:), ...
                            beads_positions(11,:), beads_positions(11,:)}, ...
                            ent_norm_abs, norm_color, true, 5, true);

                % Q3-τ1
                bezier(ax, {beads_positions(4,:), beads_positions(4,:) + [0, beads_positions(11,2), 0], ...
                            beads_positions(11,:), beads_positions(11,:)}, ...
                            ent_norm_abs, norm_color, true, 5, true);
            end

            if ThreeTangle ~= 0
                % Triple connectors Q1/Q2/Q3-τ1
                % Set P2 to the same coordinates as P3 for Q1/3-τ1 connectors
                bezier(ax, {beads_positions(2,:), beads_positions(2,:) + [0, beads_positions(11,2), 0], ...
                            beads_positions(11,:), beads_positions(11,:)}, ...
                            ThreeTangle, tangle_color, true, 5, false);

                % Just a straight line for Q2-τ1
                bezier(ax, {beads_positions(3,:), beads_positions(3,:), ...
                            beads_positions(11,:), beads_positions(11,:)}, ...
                            ThreeTangle, tangle_color, true, 5, false);

                % Q3-τ1
                bezier(ax, {beads_positions(4,:), beads_positions(4,:) + [0, beads_positions(11,2), 0], ...
                            beads_positions(11,:), beads_positions(11,:)}, ...
                            ThreeTangle, tangle_color, true, 5, false);
            end
        end
    end

    % Calculate the cylinder and cone coordinates for every axis
    nTiles = 50;
    [x_cy, y_cy, z_cy] = cylinder(cylinder_radius, 20);
    cyl_pos = {x_cy, y_cy, z_cy};
    [x_cone, y_cone, z_cone] = cylinder([cone_radius, 0]);
    cone_pos = {x_cone, y_cone, z_cone*cone_height+1};
    disk_pos = {0, 0, 1};
    % The basic label position
    label_pos = {0,0,1.5};

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
        s = surf(ax, cone_pos_temp{1}+axis_position(1), ...
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
        surf(ax, disk_pos_temp{1}+axis_position(1)+xc, ...
                 disk_pos_temp{2}+axis_position(2)+yc, ...
                 disk_pos_temp{3}+axis_position(3)+zc, ...
                 'FaceColor', axis_colors(k, :), 'EdgeColor', 'none');

        % The axis label
        if viewpoint == "top" && k == 1
            label_pos_temp = {-.4, 0, 1.5};
        elseif viewpoint == "top" && k == 2
            label_pos_temp = {0, 1.7, 0};
        end
        % Draw the label
        text(label_pos_temp{1}+axis_position(1), ...
             label_pos_temp{2}+axis_position(2), ...
             label_pos_temp{3}+axis_position(3), axis_labels{k}, ...
             'Color', axcolor(1,1,:), 'FontSize', label_fontsize, ...
             'Parent', ax,'Interpreter','Latex', ...
             'HorizontalAlignment','center');

        hold on;

    end

    hold off;

    exportgraphics(gcf, Path, 'Resolution', 600);
    close;
end

function [scale] = autoscale(beads)

    % Compute the scale factor from the maximum absolute value of the beads
    % (the id bead is excluded), rounded down to multiples of 1/4 (at least 1.0)
    max_temp = 0;
    for ii = 2:numel(beads)
        max_val = max(abs(beads{ii}), [], "all");
        if max_temp < max_val
            max_temp = max_val;
        end
    end
    scale = max(1.0, floor(1/max_temp*4)/4);
end
