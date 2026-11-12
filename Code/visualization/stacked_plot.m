% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [] = stacked_plot(ParamState, Path, options)

    %!!! stacked_plot visualizes the BEADS representations of a series of
    %    2-qubit states, which are given as a parametrized state, stacked
    %    on top of each other along the z-axis. For each state, the beads
    %    are computed (via an optimized pure state decomposition) and one
    %    figure is created and saved for every requested bead type. The
    %    labels of the beads are drawn at the topmost state and the value of
    %    the parameter next to each state. The inputs ParamState
    %    and Path are required. All other parameters are
    %    optional and are specified as name-value pairs (names are
    %    case-insensitive, can be given as string or char), e.g.
    %    stacked_plot(@(z) rho(z), "fig.png", "BeadTypes", "T",
    %                 "Parameter", "p", "ParamValues", 0:0.25:1);
    %
    %    inputs: ParamState: function handle of the parametrized state, the
    %                        density operator (4x4) for the parameter value
    %                        z is returned by ParamState(z). Alternatively,
    %                        a text expression in Dirac notation (see
    %                        parse) containing the parameter named by
    %                        Parameter, e.g., "p*dens(|Phi+>)+(1-p)*eye(4)/4"
    %            Path: path to where the final images shall be saved, the
    %                  bead type is appended to the file name of each figure
    %
    %    optional name-value inputs (default value in brackets):
    %            BeadTypes: the types of correlation function beads (list
    %                       of strings or character arrays out of "T", "E",
    %                       "C", "Sigma", "Sigma_Cl", "Sigma_D", "M",
    %                       "M_Cl", "D"), one figure is created for each
    %                       (["M","T","E","C","Sigma","Sigma_Cl","Sigma_D",
    %                       "D","M_Cl"])
    %            ParamValues: the parameter values for which the states are
    %                         computed and plotted (default: 0:0.2:1)
    %            Id: shall the id bead be visualized (boolean) (false)
    %            Parameter: the name of the parameter ("z"), given without
    %                        dollar signs (they are ignored if present). It
    %                        is plotted as the LaTeX label $Parameter$ and,
    %                        if ParamState is a text expression, it is also
    %                        the name of the parameter in the expression
    %            FontSize: font size for all labels and annotations (12)
    %            EntArcs: draw connector bars scaled by the entanglement
    %                     measure of each state (boolean) (true)
    %            EntMeasure: the entanglement measure to define the radius
    %                        of the connector arcs, "tangle" or "EoF"
    %                        (entanglement of formation) ("tangle")
    %            GridPoints: the beads resolution [theta, phi] ([200, 400])
    %            OptStructure: structure which the classical beads are
    %                          optimized for (integer): -1 = entangling
    %                          structure (the overlap with the separable
    %                          structure is minimized), 0 = no structure
    %                          optimization, 1 = separable structure (the
    %                          overlap with the separable structure is
    %                          maximized) (0)
    %            MeasurementDirection: measurement direction [x y z] of the
    %                         classical beads, used if the optimal direction
    %                         is arbitrary. If empty, the direction is
    %                         requested in the command window ([])
    %            TargetQubit: the qubit which is measured for the classical
    %                         beads (1 or 2) (1)
    %
    %    outputs: none, the figures are saved to Path

    % Check inputs, ensure correct data types and set default values
    arguments
        ParamState
        Path {mustBeTextScalar}
        options.ParamValues double {mustBeVector} = 0:0.2:1
        options.BeadTypes {mustBeText} = ["M", "T", "E", "C", "Sigma", ...
                                          "Sigma_Cl", "Sigma_D", "D", "M_Cl"]
        options.Id (1,1) logical = false
        options.Parameter {mustBeTextScalar} = "z"
        options.FontSize (1,1) double {mustBePositive} = 12
        options.EntArcs (1,1) logical = true
        options.EntMeasure {mustBeTextScalar} = "tangle"
        options.GridPoints (1,2) double = [200, 400]
        options.OptStructure (1,1) double {mustBeMember(options.OptStructure, [-1 0 1])} = 0
        options.MeasurementDirection {mustBeNumeric} = []
        options.TargetQubit (1,1) double {mustBeInteger} = 1
        options.PlotEntMeasure (1,1) logical = false
    end

    % The parameter name is used without dollar signs, the dollar signs are
    % added for the LaTeX label of the plots
    ParamName = strip(erase(string(options.Parameter), "$"));
    if ParamName ~= ""
        options.Parameter = "$" + ParamName + "$";
    else
        options.Parameter = "";
    end

    % Convert a text expression into the function handle of the parametrized
    % state, the parameter is named by the parameter name
    assert(isa(ParamState, "function_handle") || ischar(ParamState) || (isstring(ParamState) && isscalar(ParamState)), ...
           "stacked_plot:invalidParamState", ...
           "ParamState must be a function handle or a text expression.");
    if ~isa(ParamState, "function_handle")
        ParamState = parse(ParamState, "Parameter", ParamName);
    end

    % Unpack the options, the parameter values are used as row vector
    BeadTypes = string(options.BeadTypes);
    ParamValues = options.ParamValues(:).';
    Id = options.Id;
    Parameter = string(options.Parameter);
    FontSize = options.FontSize;
    EntArcs = options.EntArcs;
    EntMeasure = string(options.EntMeasure);
    GridPoints = options.GridPoints;
    OptStructure = options.OptStructure;
    MeasurementDirection = options.MeasurementDirection;
    assert(isempty(MeasurementDirection) || (numel(MeasurementDirection) == 3 && norm(MeasurementDirection) > 0), ...
           "stacked_plot:invalidDirection", ...
           "MeasurementDirection must be empty or a nonzero vector [x y z].");
    TargetQubit = options.TargetQubit;
    PlotEntMeasure = options.PlotEntMeasure;

    % Check that the bead types, the entanglement measure, and the target
    % qubit are valid
    assert(~isempty(BeadTypes) && all(ismember(BeadTypes, ["T", "E", "C", "Sigma", "Sigma_Cl", "Sigma_D", "M", "M_Cl", "D"])), ...
           "stacked_plot:invalidBeadTypes", ...
           "BeadTypes must be a non-empty list of valid bead types.");
    assert(ismember(EntMeasure, ["tangle", "EoF"]), ...
           "stacked_plot:invalidEntMeasure", ...
           "EntMeasure must be ""tangle"" or ""EoF"".");
    assert(ismember(TargetQubit, [1 2]), ...
           "stacked_plot:invalidTargetQubit", ...
           "TargetQubit must be 1 or 2.");

    % Compute the beads and the entanglement measure of every state of the
    % series. The state is decomposed into pure states with optimized
    % entanglement of formation (with padding) for the bead calculation
    % Initialize the beads and entanglement measure arrays
    Beads = cell(numel(ParamValues),1);
    EntValues = zeros(numel(ParamValues),1);
    
    for ii = 1:numel(ParamValues)
        z = ParamValues(ii);
        rho = ParamState(z);

        [~, states] = opt_ent_formation(rho, true, rank(rho)^2+1-rank(rho));
        [beads, ~] = GetBeads2Qubits(states, BeadTypes, TargetQubit, GridPoints, OptStructure, MeasurementDirection);

        Beads{ii} = beads;
        if EntMeasure ~= "EoF"
            EntValues(ii) = get_concurrence(rho)^2;
        else
            EntValues(ii) = get_entanglement_of_formation(rho);
        end
    end

    % Layout parameters
    label_fontsize = FontSize;
    tolerance = 1e-2;
    camera_angle = [47,30];

    % The basic beads positions
    beads_positions = [0,-6, 0; % Id
                       -5, 0, 0; % Q1
                       5, 0, 0; % Q2
                       0, 0, 0; % {1,2}_even
                       0, 5, 0];% {1,2}_odd

    % The target qubit which gets projected in classical mutual information
    % optimizations (identical for all states)
    targetqubit = Beads{1}{22};

    % Create one figure for each requested bead type
    for bt = 1:numel(BeadTypes)
        BeadType = BeadTypes(bt);

        % Append the bead type to the path of the figure
        path_ = regexprep(Path, '\.(png|jpg|pdf|tiff)$', '_' + BeadType + '.$1');

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

            case "Sigma"
                labels = {'Id','$Q_1$','$Q_2$','$\Sigma\{1,2\}_{even}$','$\Sigma\{1,2\}_{odd}$'};

            case "Sigma_D"
                labels = {'Id','$Q_1$','$Q_2$',append('$\Sigma_D^{(',num2str(targetqubit),')}\{1,2\}_{even}$'),append('$\Sigma_D^{',num2str(targetqubit),'}\{1,2\}_{odd}$')};

            case "Sigma_Cl"
                labels = {'Id','$Q_1$','$Q_2$',append('$\Sigma_{Cl}^{(',num2str(targetqubit),')}\{1,2\}_{even}$'),append('$\Sigma_{Cl}^{',num2str(targetqubit),'}\{1,2\}_{odd}$')};

            case "M"
                labels = {'Id','$Q_1$','$Q_2$','$M\{1,2\}_{even}$','$M\{1,2\}_{odd}$'};

            case "M_Cl"
                labels = {'Id','$Q_1$','$Q_2$',append('$M_{Cl}^{(',num2str(targetqubit),')}\{1,2\}_{even}$'),append('$M_{Cl}^{',num2str(targetqubit),'}\{1,2\}_{odd}$')};

            case "D"
                labels = {'Id','$Q_1$','$Q_2$',append('$D^{(',num2str(targetqubit),')}\{1,2\}_{even}$'),append('$D^{',num2str(targetqubit),'}\{1,2\}_{odd}$')};
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
        fig.Position = [0,0,800,600];
        ax = axes();
        axes(ax);

        % Zero array, required for color assignment
        blank = 0*X;

        set(ax, 'NextPlot', 'add');

        % Load in the beads and decide what to visualize for the different
        % correlations
        for ii = 1:numel(Beads)
            % Id
            id_bead = Beads{ii}{1};

            % Q-Beads
            bead1 = Beads{ii}{2};
            bead2 = Beads{ii}{3};

            % {1,2}_even
            bilin_even_T = Beads{ii}{4};
            bilin_even_E = Beads{ii}{5};
            bilin_even_C = Beads{ii}{6};
            bilin_even_Sigma = Beads{ii}{7};
            bilin_even_SigmaCl = Beads{ii}{8};
            bilin_even_SigmaD = Beads{ii}{9};
            bilin_even_M = Beads{ii}{10};
            bilin_even_MCl = Beads{ii}{11};
            bilin_even_ECl = Beads{ii}{20};

            % {1,2}_odd
            bilin_odd_T = Beads{ii}{12};
            bilin_odd_E = Beads{ii}{13};
            bilin_odd_C = Beads{ii}{14};
            bilin_odd_Sigma = Beads{ii}{15};
            bilin_odd_SigmaCl = Beads{ii}{16};
            bilin_odd_SigmaD = Beads{ii}{17};
            bilin_odd_M = Beads{ii}{18};
            bilin_odd_MCl = Beads{ii}{19};
            bilin_odd_ECl = Beads{ii}{21};

            switch BeadType
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

                case "Sigma"
                    beads = {id_bead, bead1, bead2, ...
                             bilin_even_Sigma, bilin_odd_Sigma};

                    beads_E = {blank, blank, blank, blank, blank};

                    beads_C = {id_bead, bead1, bead2,...
                               bilin_even_Sigma, bilin_odd_Sigma};

                case "Sigma_Cl"
                    beads = {id_bead, bead1, bead2, ...
                             bilin_even_SigmaCl, bilin_odd_SigmaCl};

                    beads_E = {blank, blank, blank, blank, blank};

                    beads_C = {id_bead, bead1, bead2,...
                               bilin_even_SigmaCl, bilin_odd_SigmaCl};

                case "Sigma_D"
                    beads = {id_bead, bead1, bead2, ...
                             bilin_even_SigmaD, bilin_odd_SigmaD};

                    beads_E = {blank, blank, blank, blank, blank};

                    beads_C = {id_bead, bead1, bead2,...
                               bilin_even_SigmaD, bilin_odd_SigmaD};

                case "M_Cl"
                    beads = {id_bead, bead1, bead2, ...
                             bilin_even_MCl, bilin_odd_MCl};

                    beads_E = {blank, blank, blank, bilin_even_ECl, bilin_odd_ECl};

                    beads_C = {id_bead, bead1, bead2,...
                               bilin_even_SigmaCl, bilin_odd_SigmaCl};

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
                               bilin_even_Sigma-bilin_even_SigmaCl, bilin_odd_Sigma-bilin_odd_SigmaCl};
            end

            % Plot the individual beads
            for bb = 1:size(beads,2)
                % The z-coordinate shift for this state
                elevation = (ii-1)*3 - (numel(Beads)-1)*3/2;

                % The beads label, draw anyway
                if ii == numel(Beads) && ((bb > 4 && any(abs(beads{bb}(:)) > tolerance)) || ...
                        (bb > 1 && bb <= 4) || (bb == 1 && Id == true))
                    % These labels are smaller
                    if bb == 4 && BeadType ~= "Sigma_Cl" && BeadType ~= "Sigma_D" && ...
                            BeadType ~= "M_Cl" && BeadType ~= "D"
                        labelshift = 0.5;
                    % These labels are larger
                    elseif bb == 4 && (BeadType == "Sigma_Cl" || BeadType == "Sigma_D" || ...
                            BeadType == "M_Cl" || BeadType == "D")
                        labelshift = 0.9;
                    % Don't need any additional shift
                    else
                        labelshift = 0;
                    end
                    % Draw the bead label
                    text(ax, beads_positions(bb,1) + labelshift, ...
                        beads_positions(bb,2)+1.5*labelshift, ...
                        beads_positions(bb,3) + (ii)*3 - (numel(Beads)-1)*3/2 - 1, ...
                        labels{bb}, 'Interpreter', 'Latex', 'FontSize', label_fontsize, ...
                        'VerticalAlignment','bottom', 'HorizontalAlignment', 'center');
                end

                % Don't draw the id bead if not specified otherwise
                % Don't draw bilinear beads if they fall below the tolerance
                if (~Id && bb == 1) || (bb > 3 && all(abs(beads{bb}(:)) < tolerance))
                    continue;
                end

                % Draw the bead surface plot
                s = surf(ax, X+beads_positions(bb,1), Y+beads_positions(bb,2), ...
                    Z+beads_positions(bb,3) + elevation, beads{bb});

                % Determine the color data and assign it to the surface plot
                set(s, 'CData', get_blended_colors(beads{bb}, beads_E{bb}, beads_C{bb}, GridPoints, false), ...
                    'EdgeColor','none');

                hold on;
            end
            % Plot the z-parameter for this state
            parval = num2str(round(ParamValues(ii),2));

            text(ax, beads_positions(3,1) + 3.5, ...
                beads_positions(3,2), ...
                beads_positions(3,3) + elevation, ...
                parval, 'Interpreter', 'Latex', ...
                'FontSize', label_fontsize,'VerticalAlignment','bottom', ...
                'HorizontalAlignment', 'center');

            if BeadType == "E"
                PlotEntMeasure = true;
            else
                PlotEntMeasure = false;
            end
            if PlotEntMeasure
                meas = "$E = ";
                if EntMeasure ~= "EoF"
                    meas = "$\tau_{12} = ";
                end
                text(ax, beads_positions(3,1) - 7.2, ...
                beads_positions(3,2), ...
                beads_positions(3,3) + elevation + 0.6, ...
                append(meas,sprintf("%.2f",EntValues(ii)), "$"), 'Interpreter', 'Latex', ...
                'FontSize', 3*label_fontsize/5,'VerticalAlignment','bottom', ...
                'HorizontalAlignment', 'center', 'Rotation', -30);
            end

            hold on;

            % Draw the connector bars with radius proportional to the
            % entanglement of formation, normalize S by dividing by log(2)
            if EntArcs
                scale = 1.;
                if EntMeasure == "EoF"
                    scale = 1/log(2);
                end
                bezier(ax, {beads_positions(2,:)+[0, 0, elevation], beads_positions(2,:)+[0, 0, elevation], ...
                    beads_positions(3,:)+[0, 0, elevation], beads_positions(3,:)+[0, 0, elevation]}, ...
                    scale*EntValues(ii), [0 0 0]);
            end
            hold on;
        end

        hold on;

        % Determine the parameter label
        parlabel = "$z$";
        if Parameter ~= ""
            parlabel = Parameter;
        end

        % Add the z-parameter label
        text(ax, beads_positions(3,1) + 3.5, ...
            beads_positions(3,2), ...
            beads_positions(3,3) + (ii)*3 - (numel(Beads)-1)*3/2 - 1, ...
            parlabel, 'Interpreter', 'Latex', ...
            'FontSize', label_fontsize,'VerticalAlignment','bottom', ...
            'HorizontalAlignment', 'center');

        % Calculate the cylinder and cone coordinates for every axis
        nTiles = 50;
        [x_cy, y_cy, z_cy] = cylinder(cylinder_radius, nTiles);
        cyl_pos = {x_cy, y_cy, z_cy};
        [x_cone, y_cone, z_cone] = cylinder([cone_radius, 0], nTiles);
        cone_pos = {x_cone, y_cone, z_cone*cone_height+1};
        disk_pos = {0, 0, 1};
        % The basic label position
        label_pos = {0,0,1.8};

        % Material properties
        shading interp;
        lightangle(100,45);
        material([1 0 0.2]);
        view(camera_angle);

        % Prevent perspective distortion
        axis equal;

        % The axis limits
        set(ax,'xlim', [-10, 6]);
        set(ax, 'ylim', [-7, 7]);
        set(ax, 'zlim', [-(numel(Beads)-1)*3/2-2, (numel(Beads)-1)*3/2+2]);

        % Hide coordinate axes
        axis off;

        axes(ax)
        % Draw the coordinate axes
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

            hold on;

            % Draw the label
            text(ax, label_pos_temp{1}+axis_position(1), ...
                label_pos_temp{2}+axis_position(2), ...
                label_pos_temp{3}+axis_position(3), axis_labels{k}, ...
                'Color', axis_colors(k, :), 'FontSize', label_fontsize, ...
                'Parent', gca,'Interpreter','Latex', ...
                'HorizontalAlignment','center');

            hold on;
        end

        % Save the figure for this bead type and close it
        exportgraphics(gcf,path_,'Resolution',600);
        close;
    end
end

function [permuted_coord] = PermuteCoordinates(coordinates, permuted_pos)

    %!!! Helper function to permute coordinates for drawing coordinate
    %    axes. The elements of the coordinate cell array are rearranged
    %    according to the given permutation
    %
    %    inputs: coordinates: the cell array of coordinates
    %            permuted_pos: the permutation, the ii-th entry gives the
    %                          index of the original coordinate which is
    %                          placed at position ii
    %
    %    outputs: permuted_coord: the cell array of the permuted coordinates

    % Initialize the coordinate array
    permuted_coord = cell(1,numel(permuted_pos));

    % Now rearrange the array according to the desired permutation
    for ii = 1 : numel(permuted_coord)
        permuted_coord{ii} = coordinates{permuted_pos(ii)};
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
