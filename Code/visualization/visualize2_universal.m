% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [] = visualize2_universal(Beads, CameraAngle, BeadType, Id, EntanglementArcs, EntMeasure, EntVal, GridPoints, RadialPlot, Normalize, Path, FontSize)

    %!!! visualize2_universal visualizes the BEADS representation of a
    %    2-qubit state for one bead type with a spherical/radial plot. The
    %    absolute value of the correlation function corresponds to the color
    %    or the radius of the beads from the origin. This function only
    %    handles 2 qubit systems
    %
    %    inputs: Beads: the beads data
    %            CameraAngle: position of the camera / viewpoint
    %            BeadType: the type of correlation function beads (possible:
    %                      "T", "E", "C", "Sigma", "Sigma_Cl", "Sigma_D",
    %                      "D", "M", "M_Cl", "E_Cl", "E_D")
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
    %
    %    outputs: none, the figure is saved to Path

    % Resolution of plot
    grid_points_theta = GridPoints(1);
    grid_points_phi = GridPoints(2);

    % Layout data
    label_shift = 0.7;
    label_fontsize = FontSize;
    tolerance = 1e-4;

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

    % Determine the beads labels
    switch BeadType
        case "T"
            labels = {'Id','$Q_1$','$Q_2$','$T\{1,2\}_{even}$','$T\{1,2\}_{odd}$'};

        case "E"
            labels = {'Id','$Q_1$','$Q_2$','$E\{1,2\}_{even}$','$E\{1,2\}_{odd}$'};

        case "C"
            labels = {'Id','$Q_1$','$Q_2$','$C\{1,2\}_{even}$','$C\{1,2\}_{odd}$'};

        case "Sigma"
            labels = {'Id','$Q_1$','$Q_2$','$\Sigma\{1,2\}_{even}$','$\Sigma\{1,2\}_{odd}$'};

        case "Sigma_Cl"
            labels = {'Id','$Q_1$','$Q_2$',append('$\Sigma_{Cl}^{(',num2str(targetqubit),')}\{1,2\}_{even}$'), ...
                      append('$\Sigma_{Cl}^{(',num2str(targetqubit),')}\{1,2\}_{odd}$')};

        case "E_Cl"
            labels = {'Id','$Q_1$','$Q_2$',append('$E_{Cl}^{(',num2str(targetqubit),')}\{1,2\}_{even}$'), ...
                      append('$E_{Cl}^{(',num2str(targetqubit),')}\{1,2\}_{odd}$')};

        case "Sigma_D"
            labels = {'Id','$Q_1$','$Q_2$',append('$\Sigma_{D}^{(',num2str(targetqubit),')}\{1,2\}_{even}$'), ...
                      append('$\Sigma_{D}^{(',num2str(targetqubit),')}\{1,2\}_{odd}$')};

        case "E_D"
            labels = {'Id','$Q_1$','$Q_2$',append('$E_{D}^{(',num2str(targetqubit),')}\{1,2\}_{even}$'), ...
                      append('$E_{D}^{(',num2str(targetqubit),')}\{1,2\}_{odd}$')};

        case "D"
            labels = {'Id','$Q_1$','$Q_2$',append('$D^{(',num2str(targetqubit),')}\{1,2\}_{even}$'), ...
                      append('$D^{(',num2str(targetqubit),')}\{1,2\}_{odd}$')};

        case "M"
            labels = {'Id','$Q_1$','$Q_2$','$M\{1,2\}_{even}$','$M\{1,2\}_{odd}$'};

        case "M_Cl"
            labels = {'Id','$Q_1$','$Q_2$',append('$M_{Cl}^{(',num2str(targetqubit),')}\{1,2\}_{even}$'), ...
                      append('$M_{Cl}^{(',num2str(targetqubit),')}\{1,2\}_{odd}$')};
    end

    % Coordinate axes layout parameters
    cylinder_radius = 0.05;
    cone_height = 0.2;
    cone_radius = 0.1;
    axis_position = [-9.5, 0, 0];

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

    % One axis object
    ax = axes();
    axes(ax);
    % Facilitate plotting into ax with external functions (bezier)
    set(ax, 'NextPlot', 'add');

    % Zero array, required for color assignment
    blank = 0*X;

    % Assign the data for plotting
    switch BeadType
        % Total correlation function
        case "T"
            beads = {id_bead, bead1, bead2, ...
                     bilin_even_T, bilin_odd_T};

            beads_E = {blank, blank, blank, ...
                       bilin_even_E, bilin_odd_E};

            beads_C = {id_bead, bead1, bead2, ...
                       bilin_even_C, bilin_odd_C};

        % Entangling structure
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

        % Part of entangling structure contributing to classical
        % correlation
        case "E_Cl"
            beads = {id_bead, bead1, bead2, ...
                     bilin_even_ECl, bilin_odd_ECl};

            beads_E = {blank, blank, blank, ...
                       bilin_even_ECl, bilin_odd_ECl};

            beads_C = {id_bead, bead1, bead2, ...
                       blank, blank};

        % Part of separable structure contributing to classical correlation
        case "Sigma_Cl"
            beads = {id_bead, bead1, bead2, ...
                     bilin_even_SigmaCl, bilin_odd_SigmaCl};

            beads_E = {blank, blank, blank, blank, blank};

            beads_C = {id_bead, bead1, bead2,...
                       bilin_even_SigmaCl, bilin_odd_SigmaCl};

        % Nonclassical part of entangling structure
        case "E_D"
            beads = {id_bead, bead1, bead2, ...
                     bilin_even_E-bilin_even_ECl, bilin_odd_E-bilin_odd_ECl};

            beads_E = {blank, blank, blank, ...
                       bilin_even_E-bilin_even_ECl, bilin_odd_E-bilin_odd_ECl};

            beads_C = {id_bead, bead1, bead2, ...
                       blank, blank};

        % Nonclassical correlation structure
        case "D"
            beads = {id_bead, bead1, bead2, ...
                     bilin_even_M-bilin_even_MCl, bilin_odd_M-bilin_odd_MCl};

            beads_E = {blank, blank, blank, ...
                       bilin_even_E-bilin_even_ECl, bilin_odd_E-bilin_odd_ECl};

            beads_C = {id_bead, bead1, bead2, ...
                       bilin_even_SigmaD, bilin_odd_SigmaD};

        % Nonclassical correlations arising from separable structure
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

    end

    % Normalize the beads, i.e., scale the radius (radial plot) or the color
    % data (otherwise) by the automatically determined scale factor
    radius_scale = 1;
    color_scale = 1;
    if Normalize
        if RadialPlot
            radius_scale = 1.0*autoscale(beads);
        else
            color_scale = autoscale(beads);
        end
    end

    if Normalize
        if radius_scale ~= 1.0 || color_scale ~= 1.0
            if RadialPlot
                annotation("textbox",[0.3, 0.95, 0.057, 0.0425], 'String', append('$\times$', sprintf('%.2f', radius_scale)), 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
            else
                annotation("textbox",[0.315, 0.73, 0.07, 0.0425], 'String', append('$\times$', sprintf('%.2f', color_scale)), 'Interpreter','latex','HorizontalAlignment','center', 'FontSize',FontSize, 'LineStyle','none');
            end
        end
    end

    % Plot the individual beads (id bead, Q-beads, and the bilinear even and
    % odd beads)
    for bb = 1:size(beads,2)
        % The beads label
        if viewpoint ~= "none" && ...
           ((bb > 3 && any(abs(beads{bb}(:)) > tolerance)) || ...
           (bb > 1 && bb <= 3) || (bb == 1 && Id && Normalize == false))
            % Draw the bead label
            text(ax, beads_positions(bb,1)+label_positions(1), ...
                     beads_positions(bb,2)+label_positions(2), ...
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
    set(ax,'xlim', [-10, 6]);
    set(ax, 'ylim', [-7, 7]);
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

    exportgraphics(gcf,Path,'Resolution',300);
    close;
end

function [scale] = autoscale(beads)

    % Compute the scale factor from the maximum absolute value of the beads
    % (the id bead is excluded), rounded down to multiples of 1/4 (at least
    % 1.0), or to two decimals if this is an exact match
    max_temp = 0;
    for ii = 2:numel(beads)
        max_val = max(abs(beads{ii}), [], "all");
        if max_temp < max_val
            max_temp = max_val;
        end
    end
    scale = max(1.0, floor(1/max_temp*4)/4);
    rounded_scale = round(1/max_temp,2);
    if abs(rounded_scale - 1/max_temp) < 1e-8
        scale = rounded_scale;
    end
end
