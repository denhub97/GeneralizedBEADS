% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [] = visualize_purestates2(PureStates, CameraAngle, BeadType, Id, Lines, Weights, EntDisplay, EntArcs, EntMeasure, GridPoints, Path, FontSize)

    %!!! visualize_purestates2 visualizes the BEADS representations of all
    %    pure states constituting a decomposition of a mixed 2-qubit state,
    %    each pure state in an individual subplot, for one bead type (T, E,
    %    or C). The coordinate axes are drawn in the first subplot. Optional
    %    separator lines, the weights of the pure states and their
    %    entanglement can be added
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
    %            EntDisplay: determines whether the entanglement entropy or
    %                        the two-tangle of each pure state is shown in
    %                        the top-left corner (boolean)
    %            EntArcs: draw connector bars scaled by the specific
    %                     tangles (boolean)
    %            EntMeasure: determines which entanglement measure is used
    %                        for scaling the connector lines and which will
    %                        be displayed as an annotation ("tangle" or
    %                        entropy)
    %            GridPoints: grid dimensions for plotting [theta, phi]
    %            Path: path to where the final image shall be saved
    %            FontSize: font size for all labels and annotations
    %
    %    outputs: none, the figure is saved to Path

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
    beads_positions = [0,-6, 0; % Id
                      -5, 0, 0; % Q1
                       5, 0, 0; % Q2
                       0, 0, 0; % {1,2}_even
                       0, 5, 0];% {1,2}_odd

    % Determine the viewpoint and label positions
    % Top view
    if CameraAngle(1) == 0 && CameraAngle(2) == 90
        viewpoint = "top";

        % Label positions
        label_positions = [label_shift,        label_shift*1.2-6,  0;
                           label_shift*1.5-5,  0,                  0;
                           label_shift*1.5+5,  0,                  0;
                           label_shift,        label_shift*1.2,    0;
                           label_shift,        label_shift*1.2+5,  0];

    % Oblique view
    elseif CameraAngle(1) == 47 && CameraAngle(2) == 30
        viewpoint = "oblique";

        % Label positions
        label_positions = [label_shift*1.5,      label_shift-6,  -0.2;
                           label_shift*1.8-4.5,  0,               0;
                           label_shift*1.8+5.5,  0,               0;
                           label_shift*1.5,      label_shift,    -0.2;
                           label_shift*1.5,      label_shift+5,  -0.2];

    % Other views we don't have labels
    else
        viewpoint = "none";
    end

    % Determine the beads labels
    switch BeadType
        % Total correlation function
        case "T"
            labels = {'Id','$Q_1$','$Q_2$','$T\{1,2\}_{even}$','$T\{1,2\}_{odd}$'};

        % Connected (entanglement-based) correlation function
        case "E"
            labels = {'Id','$Q_1$','$Q_2$','$E\{1,2\}_{even}$','$E\{1,2\}_{odd}$'};

        % Compound correlation function
        case "C"
            labels = {'Id','$Q_1$','$Q_2$','$C\{1,2\}_{even}$','$C\{1,2\}_{odd}$'};
    end

    % Coordinate axes layout parameters
    cylinder_radius = 0.05;
    cone_height = 0.2;
    cone_radius = 0.1;
    axis_position = [-5, 5, 0];
    if viewpoint == "top"
        axis_position = [-5, 3, 0];
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
    radius = 1;

    % Cartesian from spherical coordinates for 3D meshgrid
    X = cos(Phi).*sin(Theta).*radius;
    Y = sin(Phi).*sin(Theta).*radius;
    Z = cos(Theta)*radius;

    % Create the figure and subplot
    fig = figure;
    if isscalar(PureStates)
        nCols = 1;
    else
        nCols = 2;
    end
    nRows = ceil(states/nCols);
    fig.Position = [0,0,800,600*nRows/nCols];
    gap = [0. 0.];
    if viewpoint == "top"
        gap = [.02 .0];
    end
    marg_h = [.01 .01];
    marg_w = [.01 .01];
    ha = tight_subplot(nRows, nCols, gap, marg_h, marg_w);

    % Zero array, required for color assignment
    blank = 0*X;

    % Weights, entanglement entropies and two-tangles of the pure states,
    % filled in the loop over the pure states
    probabilities = zeros(states,1);
    entropies = zeros(states,1);
    two_tangles = zeros(states,1);

    % Plot each pure state in an individual subplot

    for i = 1:states
        % Norm of the weighted pure state
        norm_state = norm(PureStates{i});
        probabilities(i) = norm_state^2;
        % Renormalize pure state and calculate density operator
        rho_pure = (PureStates{i}/norm_state)*(PureStates{i}/norm_state)';
        % Collect entanglement entropies of pure states
        entropies(i) = get_von_neumann_entropy(partial_trace(rho_pure,1));
        % Collect two-tangles of pure states
        two_tangles(i) = get_concurrence(rho_pure)^2;
        % First calculate all beads for total correlation function
        % Get scaled coefficients of Lisa tensor basis operators
        scaled_coefficients = {};
        [scaled_coeff]=traverse_and_multiply(LiSAReal{2},LiSAScalings{2}, ...
                                             scaled_coefficients,rho_pure);

        % Scale spherical harmonics by corresponding coefficients
        % Id bead
        id_bead = scaled_coeff{1}*Y00;

        % Q-Beads 1 & 2
        bead1 = scaled_coeff{2}*Y1n1 + scaled_coeff{3}*Y10 + scaled_coeff{4}*Y11;
        bead2 = scaled_coeff{5}*Y1n1 + scaled_coeff{6}*Y10 + scaled_coeff{7}*Y11;

        % Bilinear terms, total correlation function
        % Scale spherical harmonics by corresponding coefficients
        % j = 0, m = 0 and j = 2, m = [-2, ..., 2] terms
        bilinear_even = scaled_coeff{8}*Y00 + scaled_coeff{12}*Y2n2 + ...
                        scaled_coeff{13}*Y2n1 + scaled_coeff{14}*Y20 + ...
                        scaled_coeff{15}*Y21 + scaled_coeff{16}*Y22;

        % j = 1, m = [-1, ..., 1] terms
        bilinear_odd = scaled_coeff{9}*Y1n1 + scaled_coeff{10}*Y10 + scaled_coeff{11}*Y11;

        % Calculate rho connected and rho compound for each pure state
        [rho_connected,rho_compound] = rho_connected_pure(rho_pure,dim);

        % Bilinear terms, connected (entanglement-based) correlation function
        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_E = {};
        [scaled_coeff_E] = traverse_and_multiply(LiSAReal{2}, LiSAScalings{2},...
                                                 scaled_coefficients_E, rho_connected);

        % Scale spherical harmonics by corresponding coefficients
        % j = 0, m = 0 and j = 2, m = [-2, ..., 2] terms
        bilinear_even_E = scaled_coeff_E{8}*Y00 + scaled_coeff_E{12}*Y2n2 + ...
                          scaled_coeff_E{13}*Y2n1 + scaled_coeff_E{14}*Y20 + ...
                          scaled_coeff_E{15}*Y21 + scaled_coeff_E{16}*Y22;

        % j = 1, m = [-1, ..., 1] terms
        bilinear_odd_E = scaled_coeff_E{9}*Y1n1 + scaled_coeff_E{10}*Y10 + scaled_coeff_E{11}*Y11;

        % Bilinear terms, compound correlation function
        % Get scaled coefficients of Lisa tensor basis
        scaled_coefficients_C = {};
        [scaled_coeff_C] = traverse_and_multiply(LiSAReal{2}, LiSAScalings{2}, ...
                                                 scaled_coefficients_C, rho_compound);

        % Scale spherical harmonics by corresponding coefficients
        % j = 0, m = 0 and j = 2, m = [-2, ..., 2] terms
        bilinear_even_C = scaled_coeff_C{8}*Y00 + scaled_coeff_C{12}*Y2n2 + ...
                          scaled_coeff_C{13}*Y2n1 + scaled_coeff_C{14}*Y20 + ...
                          scaled_coeff_C{15}*Y21 + scaled_coeff_C{16}*Y22;

        % j = 1, m = [-1, ..., 1] terms
        bilinear_odd_C = scaled_coeff_C{9}*Y1n1 + scaled_coeff_C{10}*Y10 + scaled_coeff_C{11}*Y11;

        % Determine the beads' spherical function values (beads) and values
        % used for determination of color (beads_E/C)
        switch BeadType
            % Total correlation function
            case "T"
                beads = {id_bead, bead1, bead2, ...
                         bilinear_even, bilinear_odd};

                beads_E = {blank, blank, blank, ...
                           bilinear_even_E, bilinear_odd_E};

                beads_C = {id_bead, bead1, bead2, ...
                           bilinear_even_C, bilinear_odd_C};

            % Connected (entanglement-based) correlation function
            case "E"
                beads = {id_bead, bead1, bead2, ...
                         bilinear_even_E, bilinear_odd_E};

                beads_E = {blank, blank, blank, ...
                           bilinear_even_E, bilinear_odd_E};

                beads_C = {id_bead, bead1, bead2, ...
                           blank, blank};

            % Compound correlation function
            case "C"
                beads = {id_bead, bead1, bead2, ...
                       bilinear_even_C, bilinear_odd_C};

                beads_E = {blank, blank, blank, blank, blank};

                beads_C = {id_bead, bead1, bead2,...
                           bilinear_even_C, bilinear_odd_C};
        end

        % Set the axis handle of the i-th state to subplot i
        axes(ha(i));

        % Plot the individual beads
        for bb = 1:size(beads,2)
            % Don't draw the id bead if not specified otherwise
            % Don't draw bilinear beads if they fall below the tolerance
            if (~Id && bb == 1) || (bb > 3 && all(abs(beads{bb}(:)) < tolerance))
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
            set(gca,'xlim', [-6, 6]);
            set(gca, 'ylim', [-7, 7]);
            set(gca, 'zlim', [-2, 2]);

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
        % entanglement of formation, normalize entropy by dividing by
        % log(2)
        if EntArcs
            radius = entropies(i)/log(2);
            if EntMeasure == "tangle"
                radius = two_tangles(i);
            end
            bezier(ha(i), {beads_positions(2,:), beads_positions(2,:), ...
                        beads_positions(3,:), beads_positions(3,:)}, ...
                        radius, [0 0 0]);
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
            surf(gca, cyl_pos_temp{1}+axis_position(1), ...
                      cyl_pos_temp{2}+axis_position(2), ...
                      cyl_pos_temp{3}+axis_position(3), ...
                     'FaceColor', axis_colors(k, :), 'EdgeColor', 'none');

            hold on;

            % The cone surface
            surf(gca, cone_pos_temp{1}+axis_position(1), ...
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
            surf(disk_pos_temp{1}+axis_position(1)+xc, ...
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
                 'Parent', gca,'Interpreter','Latex', ...
                 'HorizontalAlignment','center');

            hold on;
        end
    end

    hold off

    % Hide last subplot axis if density operator has odd rank
    if mod(states, 2)
        axes(ha(numel(ha)))
        axis off;
    end

    % Optionally draw lines between subplots

    if Lines
        % Draw vertical lines between columns
        if nCols > 1
            for c = 1:nCols-1
                % Get the position of a subplot in this column (first row)
                pos = get(ha(c), 'Position');
                % X position of the right edge of this subplot
                x = pos(1) + pos(3) + gap(2)/2;
                annotation(fig, 'line', [x x], [0 1], 'Color', 'k', 'LineWidth', 0.5);
            end
        end

        % Draw horizontal lines between rows
        if nRows > 1
            for r = 1:nRows-1
                % First column of this row
                idx = (r-1)*nCols + 1;
                pos = get(ha(idx), 'Position');
                y = pos(2) - gap(1)/2;
                annotation(fig, 'line', [0 1], [y y], 'Color', 'k', 'LineWidth', 0.5);
            end
        end
    end

    % If desired, add text annotations in the top-left corner of each subplot
    if Weights && numel(PureStates) > 1
        for pp = 1:states
            % Position of subplot handle
            sp_pos = get(ha(pp), 'Position');

            % Top-left corner (adjust x/y slightly to stay inside)
            x_shift = 0.075;
            y_shift = 0.075;
            if viewpoint == "top"
                x_shift = 0.05;
                y_shift = 0.06;
            end

            x = sp_pos(1) + x_shift;
            if nRows == 1
                y = 1 - y_shift;
            else
                y = (nRows - floor((pp-1)/nRows))*(1/nRows) - y_shift;
            end

            if EntDisplay && EntMeasure ~= "tangle"
                value_string = sprintf('$p_%d$ = %.2f\\%%, $S_{E}$ = %.2f', pp, probabilities(pp)*100, entropies(pp));
            elseif EntDisplay && EntMeasure == "tangle"
                value_string = sprintf('$p_%d$ = %.2f\\%%, $\\tau_{12}$ = %.2f', pp, probabilities(pp)*100, two_tangles(pp));
            else
                value_string = sprintf('$p_%d$ = %.2f\\%%', pp, probabilities(pp)*100);
            end

            % Add the probability in the top-left corner
            annotation('textbox', [x y 0.05 0.05], ...
                       'String', value_string, ...
                       'FitBoxToText', 'on', ...
                       'LineStyle', 'none', ...
                       'VerticalAlignment', 'top', ...
                       'HorizontalAlignment', 'left', ...
                       'FontSize', label_fontsize, ...
                       'Interpreter', 'latex');
        end
    end

    % Export graphic
    exportgraphics(gcf,Path,'Resolution',600);

    % Finally close the figure
    close;
end

function [PermutedCoord] = PermuteCoordinates(Coordinates, PermutedPos)

    %!!! Helper function to permute coordinates for drawing coordinate
    %    axes. The elements of the coordinate cell array are rearranged
    %    according to the given permutation
    %
    %    inputs: Coordinates: the cell array of coordinates
    %            PermutedPos: the permutation, the ii-th entry gives the
    %                         index of the original coordinate which is
    %                         placed at position ii
    %
    %    outputs: PermutedCoord: the cell array of the permuted coordinates

    % Initialize the coordinate array
    PermutedCoord = cell(1,numel(PermutedPos));

    % Now rearrange the array according to the desired permutation
    for ii = 1 : numel(PermutedCoord)
        PermutedCoord{ii} = Coordinates{PermutedPos(ii)};
    end
end

function [Xr, Yr, Zr] = get_filled_circle(V, r, N)

    %!!! Helper function to generate the surface coordinates of a filled
    %    circle with N segments and radius r perpendicular to a given
    %    vector V, centered at the origin
    %
    %    inputs: V: the vector which is perpendicular to the circle plane
    %            r: the radius of the circle
    %            N: the number of grid points in radial and angular direction
    %
    %    outputs: Xr/Yr/Zr: the N x N coordinate matrices of the rotated circle

    % Normalize the direction vector
    V = V(:);
    V = V / norm(V);

    % Create a filled circle in xy plane using polar coordinates
    theta = linspace(0, 2*pi, N);
    rho = linspace(0, r, N);
    [T, R] = meshgrid(theta, rho);

    % Just a basic circle in xy plane
    X = R .* cos(T);
    Y = R .* sin(T);
    Z = zeros(size(X));

    % Flatten grid for rotation
    pts = [X(:), Y(:), Z(:)]';

    % Compute rotation matrix from Z-axis to vector V
    z_axis = [0; 0; 1];
    if all(abs(V - z_axis) < 1e-10)
        Rmat = eye(3);
    else
        rot_axis = cross(z_axis, V);
        rot_axis = rot_axis / norm(rot_axis);
        angle = acos(dot(z_axis, V));

        % Rodrigues' rotation matrix
        K = [0,            -rot_axis(3),   rot_axis(2);
             rot_axis(3),   0,            -rot_axis(1);
            -rot_axis(2),   rot_axis(1),   0];
        Rmat = eye(3) + sin(angle)*K + (1 - cos(angle))*(K^2);
    end

    % Rotate the grid and get final coordinates
    rotated_pts = Rmat * pts;
    Xr = reshape(rotated_pts(1, :), N, N);
    Yr = reshape(rotated_pts(2, :), N, N);
    Zr = reshape(rotated_pts(3, :), N, N);
end
