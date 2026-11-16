% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [] = visualize(DensityOpExpr, Path, options)

    %!!! Main function for visualizing mixed quantum states of bi- or
    %    tripartite systems. The two inputs DensityOpExpr and Path are required.
    %    All other parameters are optional and are specified as name-value
    %    pairs (names are case-insensitive, can be given as string or char),
    %    e.g. visualize(rho, "fig.png", "PureStateVis", true, ...
    %                   "GridPoints", [100 200]);
    %   
    %    input: DensityOpExpr: The density operator describing the state,
    %                      either as a numerical matrix (4x4 or 8x8) or as
    %                      a text expression in Dirac notation (see parse)
    %           Path: The path where the figures shall be saved to
    %
    %    optional name-value inputs (default value in brackets):
    %           PureStateVis: Determines whether the pure state 
    %                      decomposition shall also be visiualized (false)
    %           BeadTypes: Array of desired bead types
    %                      (["T","E","C","M","Sigma","M_Cl","Sigma_Cl","D","Sigma_D"])
    %           CameraAngle: The camera perspective [elevation, rotation]
    %                      ([47, 30] for two qubits, [54, 45] for three qubits)
    %           GridPoints: The beads resolution for plotting [theta, phi]
    %                      ([200 400])
    %           TargetQubits: The array of targetqubits for classic.
    %                      correlation projections ([1, 1, 2, 1, 2])
    %           Id: shall the id bead be visualized (boolean) (false)
    %           EntanglementArcs: Shall the entanglment arcs be drawn
    %                      (boolean) (true)
    %           RadialPlot: Shall the beads be drawn with absolute function
    %                      value as distance from the origin (false)
    %           AutoScale: Shall the function be values be scaled 
    %                      relative to globally maximum value (false)
    %           AddRadialPlot: Additionally plot the radial variant (false)
    %           AutoScaleAddRadial: Autoscale the additional radial plot
    %                      (false)
    %           PadVectors: Number of padding vectors (rank^2+1-rank)
    %           SingleView: Plot only a single view for bipartite
    %                      systems (false)
    %           OptStructure: Structure which the M_Cl-Beads are optimized
    %                      for (integer): -1 = entangling structure (the
    %                      overlap with the separable structure is
    %                      minimized), 0 = no structure optimization, 1 =
    %                      separable structure (the overlap with the
    %                      separable structure is maximized) (0)
    %           MeasurementDirection: Measurement direction [x y z] (the
    %                      Bloch vector of the projector) of the classical
    %                      beads of two-qubit states and marginals, which is
    %                      used if the optimal direction is arbitrary and
    %                      the measurement direction can be chosen freely. If
    %                      empty, the direction is requested in the command
    %                      window when it is needed ([])
    %           OptTol: Tolerance (TolX and TolFun) of the optimization of
    %                      the pure state decomposition of three-qubit
    %                      states, ignored for two qubits (1e-6)
    %           FontSize: Font size of all labels and annotations (15)
    %           Detailed: Plot the detailed set of beads for bipartite
    %                     states (requires SingleView to be false)
    %
    %           output: None, figures are plotted directly to Path

    % Check inputs, ensure correct data types and set default values
    arguments
        DensityOpExpr
        Path {mustBeTextScalar}
        options.PureStateVis (1,1) logical = false
        options.BeadTypes {mustBeText} = ["T", "E", "C", "M", "Sigma", ...
                                          "M_Cl", "Sigma_Cl", "E_Cl", "D", "Sigma_D", "E_D"]
        options.CameraAngle (1,2) double = [NaN NaN]
        options.GridPoints (1,2) double = [200 400]
        options.TargetQubits (1,:) double = [1, 1, 2, 1, 2]
        options.Id (1,1) logical = false
        options.EntanglementArcs (1,1) logical = true
        options.RadialPlot (1,1) logical = false
        options.AutoScale (1,1) logical = false
        options.AddRadialPlot (1,1) logical = false
        options.AutoScaleAddRadial (1,1) logical = false
        options.PadVectors (1,1) double = NaN
        options.SingleView (1,1) logical = false
        options.OptStructure (1,1) double {mustBeMember(options.OptStructure, [-1 0 1])} = 0
        options.MeasurementDirection {mustBeNumeric} = []
        options.OptTol (1,1) double {mustBePositive} = 1e-6
        options.FontSize (1,1) double {mustBePositive} = 15
        options.Detailed (1,1) logical = false
    end

    % Parse and validate the density operator
    DensityOp = parse(DensityOpExpr);

    % Unpack the options
    PureStateVis = options.PureStateVis;
    BeadTypes = string(options.BeadTypes);
    CameraAngle = options.CameraAngle;
    GridPoints = options.GridPoints;
    TargetQubits = options.TargetQubits;
    Id = options.Id;
    EntanglementArcs = options.EntanglementArcs;
    RadialPlot = options.RadialPlot;
    AutoScale = options.AutoScale;
    AddRadialPlot = options.AddRadialPlot;
    AutoScaleAddRadial = options.AutoScaleAddRadial;
    PadVectors = options.PadVectors;
    if isnan(PadVectors)
        PadVectors = rank(DensityOp)^2 + 1 - rank(DensityOp);
    end
    SingleView = options.SingleView;
    OptStructure = options.OptStructure;
    OptTol = options.OptTol;
    MeasurementDirection = options.MeasurementDirection;
    assert(isempty(MeasurementDirection) || (numel(MeasurementDirection) == 3 && norm(MeasurementDirection) > 0), ...
           "visualize:invalidDirection", ...
           "MeasurementDirection must be empty or a nonzero vector [x y z].");
    FontSize = options.FontSize;
    Detailed = options.Detailed;

    % Default camera angle depends on the system size
    if any(isnan(CameraAngle))
        switch size(DensityOp,1)
            case 4
                CameraAngle = [47, 30];
            case 8
                CameraAngle = [54, 45];
        end
    end

    % Compute all beads if SingleView == false for 2-Qubit systems
    if ~SingleView && size(DensityOp,1) == 4
        BeadTypes = ["T", "E", "C", "M", "Sigma", ...
                     "M_Cl", "Sigma_Cl", "D", "Sigma_D"];
    end

    % Activate unitary reshuffling with additional zero vector padding
    if PadVectors > 0
        pad = true;
    else 
        pad = false;
    end

    % Determine which beads shall be drawn for the pure state decomposition
    % visualization
    pure_correlations = string([]);
    if PureStateVis
        % Add T-beads
        if any(ismember(BeadTypes, "T")) || ...
           any(ismember(BeadTypes, "M")) || ...
           any(ismember(BeadTypes, "D")) || ...
           any(ismember(BeadTypes, "M_Cl"))
            pure_correlations(end+1) = "T";
        end
        % Add C-beads
        if any(ismember(BeadTypes, "C")) || ...
               any(ismember(BeadTypes, "Sigma")) || ...
               any(ismember(BeadTypes, "Sigma_Cl")) || ...
               any(ismember(BeadTypes, "Sigma_D"))
            pure_correlations(end+1) = "C";
        end
        %Add E-beads
        if any(ismember(BeadTypes, "E"))
            pure_correlations(end+1) = "E";
        end
    end

   % Handle the optimizations and plotting
   switch size(DensityOp,1)
        % Bipartite state
        case 4
            % Optimize decomposition wrt EoF
            [~, states] = opt_ent_formation(DensityOp, pad, PadVectors);
            % Calculate the beads
            [beads, ~] = GetBeads2Qubits(states, BeadTypes, TargetQubits(1), GridPoints, OptStructure, MeasurementDirection);
            
            % Plot all desired bead variants, name figures dynamically
            for cc = 1:numel(BeadTypes)
                % Determine the current bead type
                correlation = BeadTypes(cc);
                
                % Plot the single view and return directly
                if SingleView == false
                    if ~Detailed
                        visualize2_complete_view(beads, CameraAngle, Id, EntanglementArcs, "tangle", ...
                                                 get_concurrence(DensityOp)^2, GridPoints, RadialPlot, AutoScale, Path, FontSize, AutoScale);
                        if RadialPlot == false && AddRadialPlot == true
                            path__ = regexprep(Path, '\.(png|jpg|pdf|tiff)$', '_' + correlation + '_radial' + '.$1');
                            visualize2_complete_view(beads, CameraAngle, Id, EntanglementArcs, "tangle", ...
                                                 get_concurrence(DensityOp)^2, GridPoints, AddRadialPlot, AutoScale, path__, FontSize, AutoScale);
                        end
                    else
                        visualize2_complete_detailed(beads, CameraAngle, Id, EntanglementArcs, "tangle", ...
                                                 get_concurrence(DensityOp)^2, GridPoints, RadialPlot, AutoScale, Path, FontSize, AutoScale);
                        if RadialPlot == false && AddRadialPlot == true
                            path__ = regexprep(Path, '\.(png|jpg|pdf|tiff)$', '_' + correlation + '_radial' + '.$1');
                            visualize2_complete_detailed(beads, CameraAngle, Id, EntanglementArcs, "tangle", ...
                                                 get_concurrence(DensityOp)^2, GridPoints, AddRadialPlot, AutoScale, path__, FontSize, AutoScale);
                        end
                    end
                    break
                else
                    % Modify the path and add current bead type
                    path_ = regexprep(Path, '\.(png|jpg|pdf|tiff)$', '_' + correlation + '.$1');
                    
                    % Plot
                    visualize2_universal(beads, CameraAngle, correlation, Id, EntanglementArcs, "tangle", ...
                                             get_concurrence(DensityOp)^2, GridPoints, RadialPlot, AutoScale, path_, FontSize);
                    
                    % Plot an additional radial plot if required
                    if RadialPlot == false && AddRadialPlot == true
                        path__ = regexprep(Path, '\.(png|jpg|pdf|tiff)$', '_' + correlation + '_radial' + '.$1');
                        visualize2_universal(beads, CameraAngle, correlation, Id, EntanglementArcs, "tangle", ...
                                         get_concurrence(DensityOp)^2, GridPoints, true, (AutoScale || AutoScaleAddRadial), path__, FontSize);
                    end
                end
            end
            
            % Plot pure state decompositions if desired
            if PureStateVis
                for cc = 1:numel(pure_correlations)
                    % Get the current pure state bead type
                    pure_corr = pure_correlations(cc);
                    % Adjust the path
                    pure_path = regexprep(Path, '\.(png|jpg|pdf|tiff)$', '_pure_decomp' + pure_corr + '.$1');
                    % Plot
                    visualize_purestates2(states, CameraAngle, pure_corr, Id, true, true, true, ...
                                         true, "tangle", GridPoints, ...
                                         pure_path, FontSize);
                end
            end
        
        % Tripartite state
        case 8
            % Load the 3-qubit Pauli basis into memory
            LoadPaulis;
            % Recursive handling of subsystems
            % Optimize all subsystems, Pauli operators are handed over
            % individually, ugly but fast in Matlab
            % Tripartite optimization
            [~, states_123] = opt_connected_tripartite(DensityOp, pad, PadVectors, OptTol, ...
                                                     pauli_1x, pauli_1y, pauli_1z, ...
                                                     pauli_2x, pauli_2y, pauli_2z, ...  
                                                     pauli_3x, pauli_3y, pauli_3z, ...
                                                     pauli_1x_2x, pauli_1x_2y, pauli_1x_2z, ...
                                                     pauli_1y_2x, pauli_1y_2y, pauli_1y_2z, ...
                                                     pauli_1z_2x, pauli_1z_2y, pauli_1z_2z, ...
                                                     pauli_1x_3x, pauli_1x_3y, pauli_1x_3z, ...
                                                     pauli_1y_3x, pauli_1y_3y, pauli_1y_3z, ...
                                                     pauli_1z_3x, pauli_1z_3y, pauli_1z_3z, ...
                                                     pauli_2x_3x, pauli_2x_3y, pauli_2x_3z, ...
                                                     pauli_2y_3x, pauli_2y_3y, pauli_2y_3z, ...
                                                     pauli_2z_3x, pauli_2z_3y, pauli_2z_3z, ...
                                                     pauli_1x_2x_3x, pauli_1x_2x_3y, pauli_1x_2x_3z, ...
                                                     pauli_1x_2y_3x, pauli_1x_2y_3y, pauli_1x_2y_3z, ...
                                                     pauli_1x_2z_3x, pauli_1x_2z_3y, pauli_1x_2z_3z, ...
                                                     pauli_1y_2x_3x, pauli_1y_2x_3y, pauli_1y_2x_3z, ...
                                                     pauli_1y_2y_3x, pauli_1y_2y_3y, pauli_1y_2y_3z, ...
                                                     pauli_1y_2z_3x, pauli_1y_2z_3y, pauli_1y_2z_3z, ...
                                                     pauli_1z_2x_3x, pauli_1z_2x_3y, pauli_1z_2x_3z, ...
                                                     pauli_1z_2y_3x, pauli_1z_2y_3y, pauli_1z_2y_3z, ...
                                                     pauli_1z_2z_3x, pauli_1z_2z_3y, pauli_1z_2z_3z);
            % Bipartite marginal optimization
            [~, states_12] = opt_ent_formation(partial_trace(DensityOp,3));
            [~, states_13] = opt_ent_formation(partial_trace(DensityOp,2));
            [~, states_23] = opt_ent_formation(partial_trace(DensityOp,1));

            % Calculate the beads
            [beads, ent_norm_123] = GetBeads3Qubits(states_123, states_12, ...
                                                    states_13, states_23, BeadTypes, ...
                                                    TargetQubits, GridPoints, true, OptStructure, MeasurementDirection);
            
            % Plot all desired bead variants, name figures dynamically
            for cc = 1:numel(BeadTypes)
                % Determine the current bead type
                correlation = BeadTypes(cc);
                % Adjust the path
                path_ = regexprep(Path, '\.(png|jpg|pdf|tiff)$', '_' + correlation + '.$1');
                % Plot
                visualize3_universal(beads, CameraAngle, correlation, Id, EntanglementArcs, ...
                                  [get_concurrence(partial_trace(DensityOp,3))^2, ...
                                   get_concurrence(partial_trace(DensityOp,2))^2, ...
                                   get_concurrence(partial_trace(DensityOp,1))^2], ...
                                   get_three_tangle_mixed(states_123), ent_norm_123, ...
                                   GridPoints, RadialPlot, AutoScale, path_, FontSize, true);

                % Plot an additional radial plot
                if RadialPlot == false && AddRadialPlot == true
                    % Adjust the path
                    path__ = regexprep(Path, '\.(png|jpg|pdf|tiff)$', '_' + correlation + '_radial' + '.$1');
                    % Plot
                    visualize3_universal(beads, CameraAngle, correlation, Id, EntanglementArcs, ...
                                  [get_concurrence(partial_trace(DensityOp,3))^2, ...
                                   get_concurrence(partial_trace(DensityOp,2))^2, ...
                                   get_concurrence(partial_trace(DensityOp,1))^2], ...
                                   get_three_tangle_mixed(states_123), ent_norm_123, ...
                                   GridPoints, true, (AutoScale || AutoScaleAddRadial), path__, FontSize, true);
                end
            end
            
            % Plot pure state decompositions if required
            if PureStateVis
                for cc = 1:numel(pure_correlations)
                    % Determine the current pure state bead type
                    pure_corr = pure_correlations(cc);
                    % Adjust the path
                    pure_path = regexprep(Path, '\.(png|jpg|pdf|tiff)$', '_pure_decomp' + pure_corr + '.$1');
                    % Plot
                    visualize_purestates3(states_123, CameraAngle, pure_corr, Id, true, ...
                                          true, true, EntanglementArcs, GridPoints, true, ...
                                          pure_path, FontSize);
                end
            end
    end
end