function [] = visualize2_complete_view(Beads, CameraAngle, Id, EntanglementArcs, EntMeasure, EntVal, GridPoints, RadialPlot, Normalize, Path, FontSize, AutoScale, ScaleFactor)

    %!!! visualize2_complete_view visualizes the BEADS representation of a
    %    2-qubit state as an overview figure with a spherical/radial plot of
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
    elseif nargin < 13
        ScaleFactor = 1.0;
    end

    % Bead types shown in the individual panels of the overview, the empty
    % entries are skipped panels
    correlations = ["T", "D", "M_Cl", "", "E", "Sigma_D", "Sigma_Cl", "M", "C", "", "", "Sigma"];

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

    % The target qubit which gets projected in classical mutual information
    % optimizations
    targetqubit = Beads{22};

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
    fig.Position = [0,0,800,600];
    % Set background color to white
    fig.Color = 'w';

    % Subplot layouts: ha0 holds the overview of the T-Beads, ha_delta the
    % difference Q-Beads, and ha the panels with the individual bilinear
    % beads, unused subplots are deleted
    ha0 = tight_subplot(2,2,[0 0],[0.02 -.135],[-0.02 -.02]);
    ha_delta = tight_subplot(3,4,[0.04 .04],[0.3 0.01],[0.01 .02]);
    ha = tight_subplot(3,4,[0.04 .04],[0.02 .3],[0.02 .02]);

    for i = 1:numel(ha0)
        if i ~= 2
            delete(ha0(i));
        end
    end

    for i = 1:numel(ha_delta)
        if i ~= 2
            delete(ha_delta(i));
        end
    end

    % Panels of the grid which are not used
    skipped = [4, 10, 11];
    for i = 1:numel(ha)
        if ismember(i, skipped)
            delete(ha(i));
        else
            set(ha(i), 'Box', 'on', ...        % draw full frame
                'XColor', 'k', ...      % black x-axis
                'YColor', 'k', ...      % black y-axis
                'LineWidth', 1.5);     % optional: thicker frame
        end
    end

    % Connecting lines and arrows between the panels, which show how the
    % beads are related to each other
    color = [0, 0, 0];
    % T->M
    annotation("line", [0.125, 0.875], [0.725, 0.725], 'LineWidth',1.5, 'Color', color);
    annotation("line", [0.125, 0.125], [0.7, 0.725], 'LineWidth',1.5, 'Color', color);
    annotation("arrow", [0.875, 0.875], [0.725, 0.46], 'LineWidth',1.5, 'Color', color, 'HeadLength', 10, 'HeadWidth', 10);
    % M->Sigma
    annotation("arrow", [0.875, 0.875], [0.22, 0.26], 'LineWidth',1.5, 'Color', color, 'HeadLength', 10, 'HeadWidth', 10);
    % M->M_Cl
    annotation("arrow", [0.86, 0.73], [0.6, 0.6], 'LineWidth',1.5, 'Color', color, 'HeadLength', 10, 'HeadWidth', 10);
    annotation("line", [0.86, 0.86], [0.6, 0.46], 'LineWidth',1.5, 'Color', color);
    % C->Sigma
    annotation("arrow", [0.23, 0.77], [0.115, 0.115], 'LineWidth',1.5, 'Color', color, 'HeadLength', 10, 'HeadWidth', 10);
    % Sigma_D/Sigma_Cl->Sigma
    annotation("arrow", [0.5, 0.77], [0.135, 0.135], 'LineWidth',1.5, 'Color', color, 'HeadLength', 10, 'HeadWidth', 10);
    annotation("line", [0.375, 0.625], [0.2, 0.2], 'LineWidth',1.5, 'Color', color);
    annotation("line", [0.5, 0.5], [0.2, 0.135], 'LineWidth',1.5, 'Color', color);
    annotation("line", [0.375, 0.375], [0.2, 0.26], 'LineWidth',1.5, 'Color', color);
    annotation("line", [0.625, 0.625], [0.2, 0.26], 'LineWidth',1.5, 'Color', color);
    % Sigma->Sigma_Cl
    annotation("arrow", [0.64, 0.64], [0.155, 0.26], 'LineWidth',1.5, 'Color', color, 'HeadLength', 10, 'HeadWidth', 10);
    annotation("line", [0.64, 0.77], [0.155, 0.155], 'LineWidth',1.5, 'Color', color);
    % M->D
    annotation("arrow", [0.375, 0.375], [0.48, 0.5], 'LineWidth',1.5, 'Color', color, 'HeadLength', 10, 'HeadWidth', 10);
    annotation("line", [0.375, 0.845], [0.48, 0.48], 'LineWidth',1.5, 'Color', color);
    annotation("line", [0.845, 0.845], [0.46, 0.48], 'LineWidth',1.5, 'Color', color);
    annotation("line", [0.625, 0.625], [0.50, 0.48], 'LineWidth',1.5, 'Color', color);

    % Descriptive text annotations
    annotation("textbox",[0.1, 0.46, 0.05, 0.0425], 'String', '$=$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.1, 0.22, 0.05, 0.0425], 'String', '$+$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.23, 0.33, 0.04, 0.0425], 'String', '$+$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.48, 0.33, 0.04, 0.0425], 'String', '$+$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.73, 0.33, 0.04, 0.0425], 'String', '$=$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.4, 0.725, 0.2, 0.0425], 'String', '$- Q_1 \cdot Q_2$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.4, 0.075, 0.2, 0.0425], 'String', '$- Q_1 \cdot Q_2$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.46, 0.15, 0.04, 0.0425], 'String', '$=$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.88, 0.22, 0.04, 0.0425], 'String', '$+E$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.7, 0.613, 0.2, 0.0425], 'String', append('Proj. $\Pi_{M_{Cl}}^{(', num2str(targetqubit), ')}$'), 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.6, 0.168, 0.2, 0.0425], 'String', append('Proj. $\Pi_{M_{Cl}}^{(', num2str(targetqubit), ')}$'), 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.69, 0.491, 0.2, 0.0425], 'String', append('$M - M_{Cl}^{(', num2str(targetqubit), ')}$'), 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');

    % Bead labels
    annotation("textbox",[0.015, 0.655, 0.04, 0.0425], 'String', '\textit{\textbf{T}}', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.52, 0.655, 0.064, 0.0425], 'String', append('\textit{\textbf{M}}$_{Cl}^{(', num2str(targetqubit), ')}$') , 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.015, 0.415, 0.04, 0.0425], 'String', '\textit{\textbf{E}}', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.015, 0.175, 0.04, 0.0425], 'String', '\textit{\textbf{C}}', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.276, 0.415, 0.04, 0.0425], 'String', append('$\mathbf{\Sigma}_{D}^{(', num2str(targetqubit), ')}$'), 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.52, 0.415, 0.055, 0.0425], 'String', append('$\mathbf{\Sigma}_{Cl}^{(', num2str(targetqubit), ')}$') , 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.77, 0.415, 0.04, 0.0425], 'String', '\textit{\textbf{M}}', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.765, 0.175, 0.04, 0.0425], 'String', '$\mathbf{\Sigma}$', 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
    annotation("textbox",[0.276, 0.655, 0.04, 0.0425], 'String', append('\textit{\textbf{D}}$^{(', num2str(targetqubit), ')}$'), 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');

    % Plot boxes
    annotation("rectangle", [0.02, 0.5, 0.21, 0.2], 'LineWidth', 1.5)
    annotation("rectangle", [0.27, 0.5, 0.21, 0.2], 'LineWidth', 1.5)
    annotation("rectangle", [0.52, 0.5, 0.21, 0.2], 'LineWidth', 1.5)
    annotation("rectangle", [0.02, 0.26, 0.21, 0.2], 'LineWidth', 1.5)
    annotation("rectangle", [0.27, 0.26, 0.21, 0.2], 'LineWidth', 1.5)
    annotation("rectangle", [0.52, 0.26, 0.21, 0.2], 'LineWidth', 1.5)
    annotation("rectangle", [0.77, 0.26, 0.21, 0.2], 'LineWidth', 1.5)
    annotation("rectangle", [0.02, 0.02, 0.21, 0.2], 'LineWidth', 1.5)
    annotation("rectangle", [0.77, 0.02, 0.21, 0.2], 'LineWidth', 1.5)

    % Zero array, required for color assignment
    blank = 0*X;

    % Scaling of the radius (radial plot) or of the color data (otherwise)
    % of the T-Beads, either automatically or by the custom scale factor
    radius_scale = 1;
    color_scale = 1;
    beads = {blank, bead1, bead2, ...
        bilin_even_T, bilin_odd_T};

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

    if AutoScale || ScaleFactor ~= 1.0
        if radius_scale ~= 1.0 || color_scale ~= 1.0
            if RadialPlot
                annotation("textbox",[0.01, 0.95, 0.057, 0.0425], 'String', append('$\times$', sprintf('%.2f', radius_scale)), 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
            else
                annotation("textbox",[0.01, 0.95, 0.07, 0.0425], 'String', append('$\times$', sprintf('%.2f', color_scale)), 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
            end
        end
    end

    % Draw the delta Q-Beads if required

    if (any(abs(delta_bead1) > tolerance, 'all') || any(abs(delta_bead2) > tolerance, 'all')) && ...
       (any(abs(bilin_even_M) > tolerance, 'all') || any(abs(bilin_odd_M) > tolerance, 'all'))

        annotation("rectangle", [0.27, 0.78, 0.21, 0.2], 'LineWidth', 1.5, 'LineStyle','--')

        ax = ha_delta(2);
        hold(ax, 'on')
        delta_beads_positions = [0, 0, 0;
                                 -3, 0, 0;
                                 3, 0, 0];
        delta_beads_labels = {"", append("$\Delta^{(", num2str(targetqubit), ")}Q_1$"), append("$\Delta^{(", num2str(targetqubit), ")}Q_2$")};
        beads = {blank, -delta_bead1, -delta_bead2};
        beads_E = {blank, blank, blank};
        beads_C = {blank, -delta_bead1, -delta_bead2};
        placeholder_pos = [0.3055 0.887;
                           0.388 0.827];

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
                annotation("ellipse", ...
                    [placeholder_pos(bb-1, 1), ...
                     placeholder_pos(bb-1, 2), ...
                     .75*0.0776, ...
                     0.0776], ...
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
        delete(ha_delta(2));
    end

    % Plot the full individual beads (T-Beads picture)
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
    set(ax,'xlim', [-6, 6.2]);
    set(ax, 'ylim', [-5.5, 7]);
    set(ax, 'zlim', [-2, 2]);

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

    % If the T-Beads do not contain bilinear correlations, show the
    % M-Beads in the overview instead and adjust the scaling accordingly
    if all(abs(beads{4}(:)) < tolerance) && all(abs(beads{5}(:)) < tolerance)
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

    if AutoScale || ScaleFactor ~= 1.0
        if radius_scale ~= 1.0 || color_scale ~= 1.0
            if RadialPlot
                annotation("textbox",[0.01, 0.725, 0.057, 0.0425], 'String', append('$\times$', sprintf('%.2f', radius_scale/(1+1/3))), 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
            else
                annotation("textbox",[0.01, 0.725, 0.07, 0.0425], 'String', append('$\times$', sprintf('%.2f', color_scale)), 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
            end
        end
    end

    % Plot the bilinear beads of each bead type in its panel of the grid
    for corr = 1:numel(ha)
        if ismember(corr, skipped)
            continue
        end
        correl = correlations(corr);
        ax = ha(corr);
        hold(ax, 'on')

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
                    an_yshift = 0.028;
                    an_xshift = 0.03;
                    if bb == 5
                        an_yshift = 0.088;
                        an_xshift = 0.12;
                    end

                    annotation("ellipse", ...
                        [0.02 + (0.21+.04)*mod(corr-1, 4) + an_xshift, ...
                        0.02 + (0.2+.04)*(2-floor((corr-1)/4)) + an_yshift, ...
                        .75*0.0776, ...
                        0.0776], ...
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
