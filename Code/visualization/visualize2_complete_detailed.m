% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [] = visualize2_complete_detailed(Beads, CameraAngle, Id, EntanglementArcs, EntMeasure, EntVal, GridPoints, RadialPlot, Normalize, Path, FontSize, AutoScale, ScaleFactor, CorrelType)

    %!!! visualize2_complete_detailed visualizes the BEADS representation of a
    %    2-qubit state as a detailed overview figure with a spherical/radial plot of
    %    the total correlation beads (including the Q-beads and the
    %    coordinate axes) and of the bilinear beads of all bead types in an
    %    arrangement which shows how the beads are related to each other,
    %    i.e., T, E, C, M, the separable structure (Sigma), the classical
    %    beads (SigmaCl, MCl), and the nonclassical beads (SigmaD, D). The
    %    absolute value of the correlation function corresponds to the color
    %    or the radius of the beads from the origin. This function only
    %    handles 2 qubit systems
    %
    %    inputs: Beads: the beads data
    %            CameraAngle: position of the camera / viewpoint
    %            Id: specify whether id bead shall be drawn (boolean)
    %            EntanglementArcs: specify whether connector arcs shall be
    %                              drawn (boolean)
    %            EntMeasure: the entanglement measure to define the radius
    %                        of connector arcs
    %            EntVal: the entanglement_measure value
    %            GridPoints: the beads resolution [theta, phi]
    %            RadialPlot: specify whether the beads shall be plotted with
    %                        absolute spherical function value as radius from
    %                        the origin
    %            Normalize: Normalize the droplets, this will not draw the
    %                       id droplet (boolean)
    %            Path: path to where the final image shall be saved
    %            FontSize: font size for all labels and annotations
    %            AutoScale: (optional, default true) scale the spherical
    %                       function values automatically relative to the
    %                       globally maximum value
    %            ScaleFactor: (optional, default 1.0) scale the spherical
    %                         function values by a custom factor
    %
    %    outputs: none, the figure is saved to Path

    if nargin < 12
        AutoScale = true;
        ScaleFactor = 1.0;
        CorrelType = "";
    elseif nargin < 13
        ScaleFactor = 1.0;
        CorrelType = "";
    elseif nargin < 14
        CorrelType = "";
    end

    % The target qubit which gets projected in classical mutual information
    % optimizations
    targetqubit = Beads{22};

    % Bead types shown in the individual panels of the overview, the empty
    % entries are skipped panels
    DrawDeltaQ = false;
    DrawT = false;
    if CorrelType == "Entanglement"
        correlations = ["","","","","","","E", "E_Cl", "E_D"];
        boxlabels = ["$\textit{\textbf{E}}$", ...
                      append("\textit{\textbf{E}}$_{Cl}^{(", num2str(targetqubit), ")}$"), ...
                      append("\textit{\textbf{E}}$_{D}^{(", num2str(targetqubit), ")}$")];
    elseif CorrelType == "Classical"
        correlations = ["","","","","","","M_Cl", "E_Cl", "Sigma_Cl"];
        boxlabels = [append("\textit{\textbf{M}}$_{Cl}^{(", num2str(targetqubit), ")}$"), ...
                     append("\textit{\textbf{E}}$_{Cl}^{(", num2str(targetqubit), ")}$"), ...
                     append("$\mathbf{\Sigma}_{Cl}^{(", num2str(targetqubit), ")}$")];
    elseif CorrelType == "Nonclassical"
        correlations = ["","","","","","","D", "E_D", "Sigma_D"];
        boxlabels = [append("\textit{\textbf{D}}$^{(", num2str(targetqubit), ")}$"), ...
                     append("\textit{\textbf{E}}$_{D}^{(", num2str(targetqubit), ")}$"), ...
                     append("$\mathbf{\Sigma}_{D}^{(", num2str(targetqubit), ")}$")];
        DrawDeltaQ = true;
    elseif CorrelType == "Separable"
        correlations = ["","","","","","","Sigma", "Sigma_Cl", "Sigma_D"];
        boxlabels = ["$\mathbf{\Sigma}$", ...
                      append("$\mathbf{\Sigma}_{Cl}^{(", num2str(targetqubit), ")}$"), ...
                      append("$\mathbf{\Sigma}_{D}^{(", num2str(targetqubit), ")}$")];
    elseif CorrelType == "TotEntSep"
        correlations = ["","","","","","","M", "E", "Sigma"];
        boxlabels = ["$\textit{\textbf{M}}$", ...
                     "$\textit{\textbf{E}}$", ...
                     "$\mathbf{\Sigma}$"];
    elseif CorrelType == "TotClasNonclas"
        correlations = ["","","","","","","M", "M_Cl", "D"];
        boxlabels = ["$\textit{\textbf{M}}$", ...
                     append("\textit{\textbf{M}}$_{Cl}^{(", num2str(targetqubit), ")}$"), ...
                     append("\textit{\textbf{D}}$^{(", num2str(targetqubit), ")}$")];
        DrawDeltaQ = true;
    elseif CorrelType == "Minimal"
        correlations = ["", "", "", "", "E_Cl","E_D","","Sigma_Cl","Sigma_D",];
        boxlabels = ["", ...
                      append("$\mathbf{\Sigma}_{Cl}^{(", num2str(targetqubit), ")}$"), ...
                      append("$\mathbf{\Sigma}_{D}^{(", num2str(targetqubit), ")}$")];
    else
        correlations = ["M", "M_Cl", "D", "E", "E_Cl", "E_D", "Sigma", "Sigma_Cl", "Sigma_D"];
        boxlabels = ["$\mathbf{\Sigma}$", ...
                      append("$\mathbf{\Sigma}_{Cl}^{(", num2str(targetqubit), ")}$"), ...
                      append("$\mathbf{\Sigma}_{D}^{(", num2str(targetqubit), ")}$")];
        DrawDeltaQ = true;
        DrawT = true;
        CorrelType = "";
    end

    % Resolution of plot
    grid_points_theta = GridPoints(1);
    grid_points_phi = GridPoints(2);

    % Layout data
    label_shift = 0.7;
    label_fontsize = FontSize;
    tolerance = 1e-2;

    % The basic beads positions
    beads_positions = [0,-6, 0; % Id
        -5, 0, 0; % Q1
        5, 0, 0; % Q2
        0, 0, 0; % {1,2}_even
        0, 5, 0];% {1,2}_odd

    % Determine viewpoint and label positions
    if CameraAngle(1) == 0 && CameraAngle(2) == 90
        %Viewpoint
        viewpoint = 'top';
    elseif CameraAngle(1) == 47 && CameraAngle(2) == 30
        % Viewpoint
        viewpoint = 'oblique';
    else
        % Viewpoint
        viewpoint = 'none';
    end

    % Label positions
    label_positions = [label_shift*1.2,        label_shift*1.2,     0];

    % Assign data to the specific beads
    % Q-Beads
    id_bead = Beads{1};
    bead1 = Beads{2};
    bead2 = Beads{3};

    % Correlation function types for even and odd case in the following order:
    % T-Bead: Total correlation function
    % E-Bead: Connected correlation function (quantum entanglement)
    % C-Bead: Compound, C = T - E
    % Sigma-Bead: Genuine correlations arising from separable structure
    % SigmaCl-Bead: Separable contribution to classical correlation.
    %          (in other words: such that information on remaining subsystem is maximal)
    % SigmaD-Bead: Nonclassical separable structure
    % M-Bead: (Total) mutual information: total correlation function where reducible
    %         components of n-partite correlation functions have been subtracted
    % MCl-Bead: Classical mutual information correlation function. M-Beads are projected
    %          using optimized projectors (same as for SigmaCl-Beads)
    % ECl-Bead: Connected correlation function projected using the optimized
    %          projectors (same as for SigmaCl-Beads)

    % {1,2}_even
    bilin_even_T = Beads{4};
    bilin_even_E = Beads{5};
    bilin_even_C = Beads{6};
    bilin_even_Sigma = Beads{7};
    bilin_even_SigmaCl = Beads{8};
    bilin_even_SigmaD = Beads{9};
    bilin_even_M = Beads{10};
    bilin_even_MCl = Beads{11};
    bilin_even_ECl = Beads{20};

    % {1,2}_odd
    bilin_odd_T = Beads{12};
    bilin_odd_E = Beads{13};
    bilin_odd_C = Beads{14};
    bilin_odd_Sigma = Beads{15};
    bilin_odd_SigmaCl = Beads{16};
    bilin_odd_SigmaD = Beads{17};
    bilin_odd_M = Beads{18};
    bilin_odd_MCl = Beads{19};
    bilin_odd_ECl = Beads{21};

    delta_bead1 = Beads{23};
    delta_bead2 = Beads{24};

    % Determine the beads labels
    labels = {'Id','$Q_1$','$Q_2$','$T\{1,2\}_{even}$','$T\{1,2\}_{odd}$'};

    % Coordinate axes layout parameters
    cylinder_radius = 0.05;
    cone_height = 0.2;
    cone_radius = 0.1;
    axis_position = [0, -5, 0];

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
    radius = 1;

    % Cartesian from spherical coordinates for 3D meshgrid
    X = cos(Phi).*sin(Theta).*radius;
    Y = sin(Phi).*sin(Theta).*radius;
    Z = cos(Theta)*radius;

    % All beads in the same plot
    fig = figure;
    Fwidth = 600;
    Fheight = 630;
    fig.Position = [0,0,Fwidth,Fheight];
    % Set background color to white
    fig.Color = 'w';

    % Subplot layouts: ha0 holds the overview of the T-Beads, ha_delta the
    % difference Q-Beads, and ha the panels with the individual bilinear
    % beads, unused subplots are deleted
    nRows = 3;
    nCols = 3;
    gapx = .03;
    gapy = Fwidth/Fheight*.03;
    margx = 0.01;
    margxr = 0.03;
    margyt = .06;
    margyb = .3;
    boxwidth = (1 - margx - margxr - 2*gapx)/nCols;
    boxheight = (1 - margyt - margyb - 2*gapy)/nRows;

    ha0 = tight_subplot(2,1,[0.3 0],[-0.029 0],[boxwidth/5 1.65*boxwidth-boxwidth/5]);
    ha_delta = tight_subplot(nRows+1,nCols,[1.3*gapy gapx],[margyb-boxheight-1.65*gapy-0.02 margyt],[margx margxr]);
    ha = tight_subplot(nRows,nCols,[gapy gapx],[margyb margyt],[margx margxr]);

    extraGapRow = 0.015;
    extraGapCol = 0.02;

    for r = 2:nRows
        for c = 1:nCols
            idx = (r-1)*nCols + c;
            pos = get(ha(idx), 'Position');
            pos(2) = pos(2) - extraGapRow;
            set(ha(idx), 'Position', pos);
        end
    end

    for r = 2:nRows
        for c = 1:nCols
            idx = (r-1)*nCols + c;
            pos = get(ha_delta(idx), 'Position');
            pos(2) = pos(2) - extraGapRow;
            set(ha_delta(idx), 'Position', pos);
        end
    end

   
    for r = 1:nRows
        for c = 2:nCols
            idx = (r-1)*nCols + c;
            pos = get(ha(idx), 'Position');
            pos(1) = pos(1) + extraGapCol;
            set(ha(idx), 'Position', pos);
        end
    end

    for i = 1:numel(ha0)
        if i ~= 2
            delete(ha0(i));
        end
    end

    if CorrelType ~= ""
        delete(ha0(1));
    end

    for i = 1:numel(ha_delta)
        if i ~= 12
            delete(ha_delta(i));
        end
    end
    

    % Connecting lines and arrows between the panels, which show how the
    % beads are related to each other

    % Descriptive text annotations
    if CorrelType == ""
        annotation("textbox",[margx, margyb + 2*boxheight + 3*gapy/2 - extraGapRow/2, boxwidth, 0], 'String', '$=$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','middle')
        annotation("textbox",[margx, margyb + boxheight + gapy/2 - extraGapRow, boxwidth, 0], 'String', '$+$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','middle');
        
        annotation("textbox",[margx + boxwidth + gapx + extraGapCol/2, margyb + 2*boxheight + 3*gapy/2 - extraGapRow/2, boxwidth, 0], 'String', '$=$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','middle')
        annotation("textbox",[margx + boxwidth + gapx + extraGapCol/2, margyb + boxheight + gapy/2 - extraGapRow, boxwidth, 0], 'String', '$+$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','middle');
        
        annotation("textbox",[margx + 2*boxwidth + 2*gapx + extraGapCol, margyb + 2*boxheight + 3*gapy/2 - extraGapRow/2, boxwidth, 0], 'String', '$=$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','middle')
        annotation("textbox",[margx + 2*boxwidth + 2*gapx + extraGapCol, margyb + boxheight + gapy/2 - extraGapRow, boxwidth, 0], 'String', '$+$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','middle');
        
        annotation("textbox",[margx + boxwidth + extraGapCol/2, margyb + 3*boxheight/2 + gapy - extraGapRow, gapx, 0], 'String', '$=$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','middle');
        annotation("textbox",[margx + boxwidth + extraGapCol/2, margyb + 5*boxheight/2 + 2*gapy, gapx, 0], 'String', '$=$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','middle');
        
        annotation("textbox",[margx + 2*boxwidth + gapx + extraGapCol, margyb + 3*boxheight/2 + gapy  - extraGapRow, gapx, 0], 'String', '$+$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','middle');
        annotation("textbox",[margx + 2*boxwidth + gapx + extraGapCol, margyb + 5*boxheight/2 + 2*gapy, gapx, 0], 'String', '$+$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','middle');
    end

    if CorrelType ~= "Minimal"
        annotation("textbox",[margx + boxwidth + extraGapCol/2, margyb + boxheight/2 - extraGapRow, gapx, 0], 'String', '$=$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','middle');
        annotation("textbox",[margx + 2*boxwidth + gapx + extraGapCol, margyb + boxheight/2 - extraGapRow, gapx, 0], 'String', '$+$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','middle');
    end

    % Bead labels
    if CorrelType == ""
        annotation("textbox",[margx, margyb + 3*boxheight + 2*gapy, margx, 0], 'String', '$\textit{\textbf{M}}$', 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','top');
        annotation("textbox",[margx, margyb + 2*boxheight + gapy - extraGapRow , margx, 0], 'String', '$\textit{\textbf{E}}$', 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','top');
        annotation("textbox",[margx + boxwidth + gapx + extraGapCol, margyb + 3*boxheight + 2*gapy, margx, 0], 'String', append('\textit{\textbf{M}}$_{Cl}^{(', num2str(targetqubit), ')}$'), 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','top');
        annotation("textbox",[margx + 2*boxwidth + 2*gapx + extraGapCol, margyb + 3*boxheight + 2*gapy, margx, 0], 'String', append('\textit{\textbf{D}}$^{(', num2str(targetqubit), ')}$'), 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','top');
    end

    if CorrelType == "" || CorrelType == "Minimal"
        annotation("textbox",[margx + boxwidth + gapx + extraGapCol, margyb + 2*boxheight + gapy - extraGapRow, margx, 0], 'String', append('\textit{\textbf{E}}$_{Cl}^{(', num2str(targetqubit), ')}$'), 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','top');
        annotation("textbox",[margx + 2*boxwidth + 2*gapx + extraGapCol, margyb + 2*boxheight + gapy - extraGapRow, margx, 0], 'String', append('\textit{\textbf{E}}$_{D}^{(', num2str(targetqubit), ')}$'), 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','top');
    end
    
    annotation("textbox",[margx, margyb + boxheight - extraGapRow, margx, 0], 'String', boxlabels(1), 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','top');
    annotation("textbox",[margx + boxwidth + gapx + extraGapCol, margyb + boxheight - extraGapRow, margx, 0], 'String', boxlabels(2), 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','top');
    annotation("textbox",[margx + 2*boxwidth + 2*gapx + extraGapCol, margyb + boxheight - extraGapRow, margx, 0], 'String', boxlabels(3), 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','top');


    % Plot boxes
    if CorrelType == ""
        annotation("rectangle", [margx, margyb + 2*boxheight + 2*gapy, boxwidth, boxheight], 'LineWidth', 1.5)
        annotation("rectangle", [margx + boxwidth + gapx + extraGapCol, margyb + 2*boxheight + 2*gapy, boxwidth, boxheight], 'LineWidth', 1.5)
        annotation("rectangle", [margx + 2*boxwidth + 2*gapx + extraGapCol, margyb + 2*boxheight + 2*gapy, boxwidth, boxheight], 'LineWidth', 1.5)

        annotation("rectangle", [margx, margyb + boxheight + gapy - extraGapRow, boxwidth, boxheight], 'LineWidth', 1.5)
    end

    if CorrelType == "" || CorrelType == "Minimal"
        annotation("rectangle", [margx + boxwidth + gapx + extraGapCol, margyb + boxheight + gapy - extraGapRow, boxwidth, boxheight], 'LineWidth', 1.5)
        annotation("rectangle", [margx + 2*boxwidth + 2*gapx + extraGapCol, margyb + boxheight + gapy - extraGapRow, boxwidth, boxheight], 'LineWidth', 1.5)
    end

    if CorrelType ~= "Minimal"
        annotation("rectangle", [margx, margyb - extraGapRow, boxwidth, boxheight], 'LineWidth', 1.5)
    end

    annotation("rectangle", [margx + boxwidth + gapx + extraGapCol, margyb - extraGapRow, boxwidth, boxheight], 'LineWidth', 1.5)
    annotation("rectangle", [margx + 2*boxwidth + 2*gapx + extraGapCol, margyb - extraGapRow, boxwidth, boxheight], 'LineWidth', 1.5)

    if DrawT
        annotation("rectangle", [margx, 0.01, 2*boxwidth + gapx + extraGapCol, boxheight+2*gapy-extraGapRow], 'LineWidth', 1.5)
    end

    % Zero array, required for color assignment
    blank = 0*X;

    % Scaling of the radius (radial plot) or of the color data (otherwise)
    % of the T-Beads, either automatically or by the custom scale factor
    radius_scale = 1;
    color_scale = 1;
    if DrawT
        beads = {blank, bead1, bead2, ...
            bilin_even_T, bilin_odd_T};
    else
        beads = {blank, bead1, bead2, ...
            blank, blank};
    end

    if AutoScale
        if RadialPlot
            radius_scale = 1.0*autoscale(beads);
        else
            color_scale = autoscale(beads);
        end
    elseif ScaleFactor ~= 1.0 && AutoScale == false
        if RadialPlot
            radius_scale = ScaleFactor;
        else
            color_scale = ScaleFactor;
        end
    end

    if all(abs(beads{2}(:)) < tolerance) && ...
       all(abs(beads{3}(:)) < tolerance) && ...
       all(abs(beads{4}(:)) < tolerance) && ...
       all(abs(beads{5}(:)) < tolerance)
        radius_scale = 1.0;
        color_scale = 1.0;
    end

    if (AutoScale || ScaleFactor ~= 1.0) && (DrawDeltaQ || DrawT)
        if radius_scale ~= 1.0 || color_scale ~= 1.0
            xposition = margx-0.01;
            if ~DrawT
                xposition = margx - 0.01 + 2*gapx + extraGapCol + 2*boxwidth;
            end
            if RadialPlot
                annotation("textbox",[xposition, boxheight+2*gapy-extraGapRow + 0.003, margx, 0], 'String', append('$\times$', sprintf('%.2f', radius_scale)), 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none','VerticalAlignment','bottom');
            else
                annotation("textbox",[xposition, boxheight+2*gapy-extraGapRow + 0.003, margx, 0], 'String', append('$\times$', sprintf('%.2f', color_scale)), 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none','VerticalAlignment','bottom');
            end
        end
    end

    % Draw the delta Q-Beads if required
    if DrawDeltaQ
        if (any(abs(delta_bead1) > tolerance, 'all') || any(abs(delta_bead2) > tolerance, 'all')) && ...
           (any(abs(bilin_even_M) > tolerance, 'all') || any(abs(bilin_odd_M) > tolerance, 'all'))
    
            annotation("rectangle", [margx + 2*boxwidth + 2*gapx + extraGapCol, 0.01, boxwidth, boxheight+2*gapy-extraGapRow], 'LineWidth', 1.5, 'LineStyle','--')
    
            ax = ha_delta(12);
            hold(ax, 'on')
            delta_beads_positions = [0, 0, 0;
                                     -3, 0, 0;
                                     3, 0, 0];
            delta_beads_labels = {"", append("$\Delta^{(", num2str(targetqubit), ")}Q_1$"), append("$\Delta^{(", num2str(targetqubit), ")}Q_2$")};
            beads = {blank, -delta_bead1, -delta_bead2};
            beads_E = {blank, blank, blank};
            beads_C = {blank, -delta_bead1, -delta_bead2};
    
            for bb = 1:size(beads,2)
                if all(abs(beads{bb}) < tolerance, 'all') && bb < 2
                    continue
                end
    
                % Draw the bead label
                text(ax, delta_beads_positions(bb,1), ...
                    delta_beads_positions(bb,2)-2, ...
                    delta_beads_positions(bb,3)-2.6, ...
                    delta_beads_labels{bb}, 'Interpreter', 'Latex', 'FontSize', label_fontsize, ...
                    'VerticalAlignment','bottom', 'HorizontalAlignment', 'left');
    
                if all(abs(beads{bb}(:)) < tolerance)
                    if bb == 2
                        an_xshift = 1.04*boxwidth/5;
                        an_yshift = 5.45*boxheight/10;
                    else
                        an_xshift = 2.94*boxwidth/5;
                        an_yshift = 2.48*boxheight/10;
                    end
                    annotation("ellipse", ...
                        [margx + (2*boxwidth + 2*gapx) + extraGapCol + an_xshift - 0.02, ...
                         0.01 + an_yshift + 0.02, ...
                         0.082, ...
                         Fwidth/Fheight*0.08], ...
                        'LineStyle', ':', ...
                        'LineWidth',1.5, ...
                        'Color', [0.5, 0.5, 0.5]);
                    continue
                end
    
                % The radii array
                if RadialPlot
                    r = (1+1/2)*abs(beads{bb});
                else
                    r = 1.5*ones(grid_points_theta, grid_points_phi);
                end
    
                % Draw the bead surface plot
                if RadialPlot
                    s = surf(ax, radius_scale*r.*X+delta_beads_positions(bb,1), radius_scale*r.*Y+delta_beads_positions(bb,2), ...
                        radius_scale*r.*Z+delta_beads_positions(bb,3), beads{bb});
                else
                    s = surf(ax, r.*X+delta_beads_positions(bb,1), r.*Y+delta_beads_positions(bb,2), ...
                        r.*Z+delta_beads_positions(bb,3), radius_scale*beads{bb});
                end
    
                % Determine the color data and assign it to the surface plot
                set(s, 'CData', get_blended_colors(color_scale*beads{bb}, color_scale*beads_E{bb}, color_scale*beads_C{bb}, GridPoints, RadialPlot), ...
                    'EdgeColor','none');
    
                hold on;
            end
    
            % Material properties
            shading(ax, 'interp');
            lightangle(ax,100,45);
            if RadialPlot
                material (ax, [0.7 0.6 0.5]);
            else
                material(ax, [1 0 0.2]);
            end
            view(ax,CameraAngle);
            axis(ax, 'off')
            axis(ax, 'equal')
            % The axis limits
            set(ax,'ylim', [-2.5 1.5]);
            set(ax, 'xlim', [-4.5, 4.5]);
            if RadialPlot == false
                set(ax, 'zlim', [-1.8, 1.5]);
            else
                set(ax, 'zlim', [-1.8, 1.8]);
            end
        else
            delete(ha_delta(12));
        end
    else 
        delete(ha_delta(12));
    end

    % Plot the full individual beads (T-Beads picture)
    if DrawT
        beads = {id_bead, bead1, bead2, ...
            bilin_even_T, bilin_odd_T};
    
        beads_E = {blank, blank, blank, ...
            bilin_even_E, bilin_odd_E};
    
        beads_C = {id_bead, bead1, bead2, ...
            bilin_even_C, bilin_odd_C};
    
        ax = ha0(2);
        hold(ax, 'on')
    
        for bb = 1:size(beads,2)
            if bb == 1
                % Material properties
                shading(ax, 'interp');
                lightangle(ax,100,45);
            end
            % The beads label
            if viewpoint ~= "none" && ...
                    ((bb > 3 && any(abs(beads{bb}(:)) > tolerance)) || ...
                    (bb > 1 && bb <= 3) || (bb == 1 && Id && Normalize == false))
                % Draw the bead label
                text(ax, beads_positions(bb,1)+label_positions(1)+0.7, ...
                    beads_positions(bb,2)+label_positions(2)-0.7, ...
                    beads_positions(bb,3)+label_positions(3), ...
                    labels{bb}, 'Interpreter', 'Latex', 'FontSize', label_fontsize, ...
                    'VerticalAlignment','bottom', 'HorizontalAlignment', 'left');
            end
    
            % Don't draw the id bead if not specified otherwise
            % Don't draw bilinear beads if they fall below the tolerance
            if ((~Id || Normalize) && bb == 1) || ...
                    (bb > 3 && all(abs(beads{bb}(:)) < tolerance))
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
            set(s, 'CData', get_blended_colors(color_scale*beads{bb}, color_scale*beads_E{bb}, color_scale*beads_C{bb}, GridPoints, RadialPlot), ...
                'EdgeColor','none');
    
            hold on;
        end
    
        if RadialPlot
            material(ax, [0.7 0.6 0.5]);
        else
            material(ax, [1 0 0.2]);
        end
        view(ax,CameraAngle);
        axis(ax, 'off')
        axis(ax, 'equal')
        % The axis limits
        set(ax,'xlim', [-6, 5]);
        set(ax, 'ylim', [-3, 3]);
        set(ax, 'zlim', [0, 1]);
        set(ax, "Clipping", "off");
    
        hold on;
    
        % Draw the connector bars with radius proportional to the
        % entanglement of formation, normalize entropy by dividing by
        % log(2)
        if EntanglementArcs && EntVal >= tolerance
            if EntMeasure == "tangle"
                radius = EntVal;
            else
                radius = EntVal/log(2);
            end
            bezier(ax, {beads_positions(2,:), beads_positions(2,:), ...
                beads_positions(3,:), beads_positions(3,:)}, ...
                radius, [0 0 0]);
        end
    
        hold on;
    
        % Calculate the cylinder and cone coordinates for every axis
        nTiles = 50;
        [x_cy, y_cy, z_cy] = cylinder(cylinder_radius, nTiles);
        cyl_pos = {x_cy, y_cy, z_cy};
        [x_cone, y_cone, z_cone] = cylinder([cone_radius, 0], nTiles);
        cone_pos = {x_cone, y_cone, z_cone*cone_height+1};
        disk_pos = {0, 0, 1};
        % The basic label position
        label_pos = {0, 0, 1.5};
    
        for k = 1:3
            % Get the coordinates of this axis
            cyl_pos_temp = PermuteCoordinates(cyl_pos, axis_perm{k});
            cone_pos_temp = PermuteCoordinates(cone_pos, axis_perm{k});
            label_pos_temp = PermuteCoordinates(label_pos, axis_perm{k});
            disk_pos_temp = PermuteCoordinates(disk_pos, axis_perm{k});
    
            % The cylinder surface
            surf(ax, cyl_pos_temp{1}+axis_position(1), ...
                cyl_pos_temp{2}+axis_position(2), ...
                cyl_pos_temp{3}+axis_position(3), ...
                'FaceColor', axis_colors(k, :), 'EdgeColor', 'none');
    
            hold on;
    
            % The cone surface
            surf(ax, cone_pos_temp{1}+axis_position(1), ...
                cone_pos_temp{2}+axis_position(2), ...
                cone_pos_temp{3}+axis_position(3), ...
                'FaceColor', axis_colors(k, :), 'EdgeColor', 'none');
    
            hold on;
    
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
                'Color', axis_colors(k, :), 'FontSize', label_fontsize, ...
                'Parent', ax,'Interpreter','Latex', ...
                'HorizontalAlignment','center');
    
            hold on;
        end  
    
        hold off
    else
        delete(ha0(1));
        delete(ha0(2));
    end

    % Scale correlation function beads with respect to main bead of shown
    % data
    if CorrelType == "Entanglement"
        beads = {id_bead, bead1, bead2, ...
        bilin_even_E, bilin_odd_E};
    elseif CorrelType == "Classical"
        beads = {id_bead, bead1, bead2, ...
        bilin_even_MCl, bilin_odd_MCl};
    elseif CorrelType == "Nonclassical"
        beads = {id_bead, bead1, bead2, ...
        bilin_even_M-bilin_even_MCl, bilin_odd_M-bilin_odd_MCl};
    elseif CorrelType == "Separable"
        beads = {id_bead, bead1, bead2, ...
        bilin_even_Sigma, bilin_odd_Sigma};
    elseif CorrelType == "TotEntSep"
        beads = {id_bead, bead1, bead2, ...
        bilin_even_M, bilin_odd_M};
    elseif CorrelType == "TotClasNonclas"
        beads = {id_bead, bead1, bead2, ...
        bilin_even_M, bilin_odd_M};
    elseif CorrelType == "Minimal"
        bilin_even_ED = bilin_even_E-bilin_even_ECl;
        bilin_odd_ED = bilin_odd_E-bilin_odd_ECl;

        [~, idx] = max([max(abs(bilin_even_ECl(:))), max(abs(bilin_even_SigmaCl(:))), ...
                   max(abs(bilin_even_ED(:))), max(abs(bilin_even_SigmaD(:)))]);
        
        arrays = {bilin_even_ECl, bilin_even_SigmaCl, bilin_even_ED, bilin_even_SigmaD};
        arrays_odd = {bilin_odd_ECl, bilin_odd_SigmaCl, bilin_odd_ED, bilin_odd_SigmaD};

        beads = {id_bead, bead1, bead2, ...
                 arrays{idx}, arrays_odd{idx}};
    else
        beads = {id_bead, bead1, bead2, ...
        bilin_even_M, bilin_odd_M};
    end

    if all(abs(beads{4}(:)) < tolerance) && all(abs(beads{5}(:)) < tolerance)
        AutoScale = false;
        ScaleFactor = 1.0;
    end

    if AutoScale
        if RadialPlot
            radius_scale = (1+1/3)*autoscale(beads(1,4:end));
        else
            color_scale = autoscale(beads(1,4:end));
        end
    elseif ScaleFactor ~= 1.0 && AutoScale == false
        if RadialPlot
            radius_scale = ScaleFactor;
        else
            color_scale = ScaleFactor;
        end
    end

    if all(abs(beads{2}(:)) < tolerance) && ...
       all(abs(beads{3}(:)) < tolerance) && ...
       all(abs(beads{4}(:)) < tolerance) && ...
       all(abs(beads{5}(:)) < tolerance)
        radius_scale = 1.0;
        color_scale = 1.0;
    end

    if AutoScale || ScaleFactor ~= 1.0
        if radius_scale ~= 1.0 || color_scale ~= 1.0
            yposition = margyb - extraGapRow + 2.5*gapy + 3*boxheight - 0.007;
            xposition = margx-0.01;
            if CorrelType ~= "" && CorrelType ~= "Minimal"
                yposition = margyb - extraGapRow + boxheight - 0.007;
            elseif CorrelType == "Minimal"
                yposition = margyb - extraGapRow + 2*boxheight + gapy - 0.007;
                xposition = margx-0.01 + gapx + extraGapCol + boxwidth;
            end
            if RadialPlot
                annotation("textbox",[xposition, yposition, margx, 0], 'String', append('$\times$', sprintf('%.2f', radius_scale/(1+1/3))), 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','bottom');
            else
                annotation("textbox",[xposition, yposition, margx, 0], 'String', append('$\times$', sprintf('%.2f', color_scale)), 'Interpreter','latex','HorizontalAlignment','left', 'FontSize',FontSize, 'LineStyle','none', 'VerticalAlignment','bottom');
            end
        end
    end

    % Plot the bilinear beads of each bead type in its panel of the grid
    for corr = 1:numel(ha)
        correl = correlations(corr);
        ax = ha(corr);
        hold(ax, 'on')

        skip = false;

        % Assign the data for plotting
        switch correl
            % Total correlation function
            case "T"
                beads = {id_bead, bead1, bead2, ...
                    bilin_even_T, bilin_odd_T};

                beads_E = {blank, blank, blank, ...
                    bilin_even_E, bilin_odd_E};

                beads_C = {id_bead, bead1, bead2, ...
                    bilin_even_C, bilin_odd_C};

                % Connected (entanglement-based) correlation function
            case "E"
                beads = {id_bead, bead1, bead2, ...
                    bilin_even_E, bilin_odd_E};

                beads_E = {blank, blank, blank, ...
                    bilin_even_E, bilin_odd_E};

                beads_C = {id_bead, bead1, bead2, ...
                    blank, blank};

                % Compound correlation function
            case "C"
                beads = {id_bead, bead1, bead2, ...
                    bilin_even_C, bilin_odd_C};

                beads_E = {blank, blank, blank, blank, blank};

                beads_C = {id_bead, bead1, bead2,...
                    bilin_even_C, bilin_odd_C};

                % Correlation functions of separable components of the state
                % (classical and nonclassical correlations)
            case "Sigma"
                beads = {id_bead, bead1, bead2, ...
                    bilin_even_Sigma, bilin_odd_Sigma};

                beads_E = {blank, blank, blank, blank, blank};

                beads_C = {id_bead, bead1, bead2,...
                    bilin_even_Sigma, bilin_odd_Sigma};

                % Correlation functions of separable components of the state
                % contributing to classical correlation (CMI)
            case "Sigma_Cl"
                beads = {id_bead, bead1, bead2, ...
                    bilin_even_SigmaCl, bilin_odd_SigmaCl};

                beads_E = {blank, blank, blank, blank, blank};

                beads_C = {id_bead, bead1, bead2,...
                    bilin_even_SigmaCl, bilin_odd_SigmaCl};

                % Nonclassical correlation functions of separable components of the
                % state
            case "Sigma_D"
                beads = {id_bead, bead1, bead2, ...
                    bilin_even_SigmaD, bilin_odd_SigmaD};

                beads_E = {blank, blank, blank, blank, blank};

                beads_C = {id_bead, bead1, bead2,...
                    bilin_even_SigmaD, bilin_odd_SigmaD};

            case "E_Cl"
                beads = {id_bead, bead1, bead2, ...
                    bilin_even_ECl, bilin_odd_ECl};

                beads_E = {blank, blank, blank, bilin_even_ECl, bilin_odd_ECl};

                beads_C = {id_bead, bead1, bead2,...
                    blank, blank};

                % Nonclassical correlation functions of separable components of the
                % state
            case "E_D"
                beads = {id_bead, bead1, bead2, ...
                    bilin_even_E-bilin_even_ECl, bilin_odd_E-bilin_odd_ECl};

                beads_E = {blank, blank, blank, bilin_even_E-bilin_even_ECl, bilin_odd_E-bilin_odd_ECl};

                beads_C = {id_bead, bead1, bead2,...
                    blank, blank};

                % Classical mutual information function
            case "M_Cl"
                beads = {id_bead, bead1, bead2, ...
                    bilin_even_MCl, bilin_odd_MCl};

                beads_E = {blank, blank, blank, ...
                    bilin_even_ECl, bilin_odd_ECl};

                beads_C = {id_bead, bead1, bead2,...
                    bilin_even_SigmaCl, bilin_odd_SigmaCl};

                % Total mutual information correlation function
            case "M"
                beads = {id_bead, bead1, bead2, ...
                    bilin_even_M, bilin_odd_M};

                beads_E = {blank, blank, blank, ...
                    bilin_even_E, bilin_odd_E};

                beads_C = {id_bead, bead1, bead2, ...
                    bilin_even_Sigma, bilin_odd_Sigma};

            case "D"
                beads = {id_bead, bead1, bead2, ...
                    bilin_even_M-bilin_even_MCl, bilin_odd_M-bilin_odd_MCl};

                beads_E = {blank, blank, blank, ...
                    bilin_even_E-bilin_even_ECl, bilin_odd_E-bilin_odd_ECl};

                beads_C = {id_bead, bead1, bead2, ...
                    bilin_even_SigmaD, bilin_odd_SigmaD};

            otherwise
                skip = true;
        end

        if skip
            axis(ax, 'off')
            continue
        end

        % Plot the individual beads
        for bb = 1:size(beads,2)
            if bb < 4
                continue
            end

            % Don't draw the id bead if not specified otherwise
            % Don't draw bilinear beads if they fall below the tolerance
            if ((~Id || Normalize) && bb == 1) || ...
                    (bb > 3 && all(abs(beads{bb}(:)) < tolerance))

                if all(abs(beads{bb}(:)) < tolerance)
                    an_yshift = 0.034;
                    an_xshift = .82*boxwidth/5;
                    if bb == 5
                        an_yshift = 4.45*boxheight/10;
                        an_xshift = 2.8*boxwidth/5;
                    end
                    
                    placeholder_shift_x = 0;
                    if mod(corr,3) ~= 1
                        placeholder_shift_x = extraGapCol;
                    end
                    placeholder_shift_y = 0;
                    if corr > 3 
                        placeholder_shift_y = -extraGapRow;
                    end
                    diameter = 0.082;
                    annotation("ellipse", ...
                        [margx + (boxwidth + gapx) * mod(corr-1, 3) + an_xshift + placeholder_shift_x, ...
                         margyb + (boxheight + gapy) * (2-floor((corr-1)/3)) + an_yshift + placeholder_shift_y, ...
                        diameter, ...
                        Fwidth/Fheight*diameter], ...
                        'LineStyle', ':', ...
                        'LineWidth',1.5, ...
                        'Color', [0.5, 0.5, 0.5]);
                end
                continue;
            end

            % The radii array
            if RadialPlot
                r = abs(beads{bb});
            else
                r = 1.3*ones(grid_points_theta, grid_points_phi);
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
            set(s, 'CData', get_blended_colors(color_scale*beads{bb}, color_scale*beads_E{bb}, color_scale*beads_C{bb}, GridPoints, RadialPlot), ...
                'EdgeColor','none');

            hold on;
        end

        hold on;

        % Material properties
        shading(ax, 'interp');
        lightangle(ax,100,45);
        if RadialPlot
            material (ax, [0.7 0.6 0.5]);
        else
            material(ax, [1 0 0.2]);
        end
        view(ax,CameraAngle);
        axis(ax, 'off')
        axis(ax, 'equal')
        % The axis limits
        set(ax,'xlim', [-1.3 1.5]);
        set(ax, 'ylim', [-1.5, 6.3]);
        if RadialPlot == false
            set(ax, 'zlim', [-1.5, 1.5]);
        else
            set(ax, 'zlim', [-1.8, 1.8]);
        end

        hold off;
    end
    
    % Pixels to cut from the top
    cropTop = 15;
    % Pixels to cut from the bottom
    cropBottom = 0;     
    if CorrelType ~= ""
        if CorrelType ~= "Minimal"
            cropTop = 310;
        else
            cropTop = 170;
        end
        if ~DrawDeltaQ
            cropBottom = 175;
        end
    end

    % Old and new figure height in pixels
    fig.Units = 'pixels';
    H0 = fig.Position(4);
    H1 = H0 - cropTop - cropBottom;

    % Freeze the figure children (axes, colorbars, legends) in pixel units.
    % Only direct children are used so that text objects inside the axes
    % (data units) are not touched
    kids = fig.Children;
    kids = kids(arrayfun(@(k) isprop(k, 'Units') && isprop(k, 'Position'), kids));
    set(kids, 'Units', 'pixels');

    % Annotations (textboxes, rectangles, ellipses, arrows, lines) live in a
    % hidden overlay layer in normalized figure units and are converted
    % manually: pixel sizes stay fixed, the y-offset is shifted by cropBottom
    overlay = findall(fig, 'Type', 'annotationpane');
    ann = [];
    for o = reshape(overlay,1,[])
        ann = [ann; o.Children(:)]; %#ok<AGROW>
    end
    conv_y = @(y) (y*H0 - cropBottom)/H1;
    conv_h = @(h) h*H0/H1;
    for a = reshape(ann,1,[])
        if isprop(a, 'Units') && ~strcmp(a.Units, 'normalized')
            a.Units = 'normalized';
        end
        % Arrows and lines are defined by their end points
        if isprop(a, 'Y')
            a.Y = conv_y(a.Y);
        % Textboxes, rectangles and ellipses are defined by [x y w h]
        elseif isprop(a, 'Position')
            p = a.Position;
            a.Position = [p(1), conv_y(p(2)), p(3), conv_h(p(4))];
        end
    end

    % Shrink the figure
    fig.Position(4) = H1;

    % Pixel origin is bottom-left, so shift the axes down if the bottom is cropped
    if cropBottom > 0
        for k = reshape(kids,1,[])
            k.Position(2) = k.Position(2) - cropBottom;
        end
    end
    % Save the figure to the specified path and close it
    exportgraphics(gcf,Path,'Resolution',300);
    close;
end

function [scale] = autoscale(beads)

    % Compute the scale factor from the maximum absolute value of the beads,
    % rounded down to multiples of 1/4 (at least 1.0)
    max_temp = 0;
    for ii = 1:numel(beads)
        max_val = max(abs(beads{ii}), [], "all");
        if max_temp < max_val
            max_temp = max_val;
        end
    end
    scale = max(1.0, floor(1/max_temp*4)/4);
end
