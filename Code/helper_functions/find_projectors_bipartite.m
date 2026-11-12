% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [RhoSigmaCl, RhoMCl, RhoDCl, RhoECl, DeltaRho1, DeltaRho2] = find_projectors_bipartite(DensityOp, RhoSigma, RhoM, RhoConn, TargetQubit, OptStructure, MeasurementDirection)

    %!!! Find the single-qubit measurement projectors on the target qubit of
    %    a 2-qubit state and apply them to the correlation structures of the
    %    state. The projectors are optimized for the density operator
    %    (maximal classical correlations). If the optimal orientation of
    %    the projectors is not unique, i.e., if the optimal measurement
    %    direction has a degree of freedom (a plane or arbitrary
    %    rotations), the projectors are chosen among the optimal ones such
    %    that they best fit the separable structure of the state. The
    %    degree of freedom is determined from repeated optimizations
    %    (no unique orientation, an optimal plane of orientations, or
    %    arbitrary orientations)
    %
    %    input: DensityOp: The density operator describing the state (4x4)
    %           RhoSigma: The separable structure of the density operator
    %           RhoM: The total correlation structure of the density operator
    %           RhoConn: The connected (entanglement) structure of the
    %                    density operator (an empty matrix can be passed,
    %                    it is then replaced by a zero matrix)
    %           TargetQubit: The qubit which is measured (1 or 2)
    %           OptStructure: Structure which the projectors are optimized
    %                       for (integer): -1 = entangling structure (the
    %                       overlap with the separable structure is
    %                       minimized), 0 = no structure optimization (the
    %                       first optimum is used), 1 = separable structure
    %                       (the overlap with the separable structure is
    %                       maximized). For -1 and 1 the degree of freedom
    %                       of the optimal projector orientation is used
    %                       (default 0)
    %           MeasurementDirection: Measurement direction [x y z] (a
    %                       vector with three elements, normalized
    %                       internally) which is used if the optimal
    %                       orientation is arbitrary. If empty (default),
    %                       the direction is requested in the command
    %                       window, where an empty input selects the z-axis
    %
    %    output: RhoSigmaCl: The classical structure of the separable
    %                        structure, i.e., RhoSigma after the measurement
    %            RhoMCl: The mutual structure after the measurement
    %            RhoDCl: The part of the separable structure which is
    %                    destroyed by the measurement, RhoSigma - RhoSigmaCl
    %            RhoECl: The connected structure after the measurement
    %            DeltaRho1/DeltaRho2: The reduced density operators of qubit
    %                                 1 and 2 after the measurement of the
    %                                 density operator

    % Handle default parameters
    if nargin < 6
        OptStructure = 0;
    end
    if nargin < 7
        MeasurementDirection = [];
    end
    assert(isempty(MeasurementDirection) || (isnumeric(MeasurementDirection) && numel(MeasurementDirection) == 3 && norm(MeasurementDirection) > 0), ...
           "find_projectors_bipartite:invalidDirection", ...
           "MeasurementDirection must be empty or a nonzero vector with three elements.");
    if isempty(RhoConn)
        RhoConn = zeros(size(DensityOp));
    end

    % Define Paulis
    pauli_x = [0 1; 1 0];
    pauli_y = [0 -1i; 1i 0];
    pauli_z = [1 0; 0 -1];

    tolerance = 1e-6;

    % Use DensityOp to get first pair of P0 and P1
    [cost,P0,P1] = opt_classcorr_bipartite(DensityOp,TargetQubit);

    % Determine the degree of freedom of the optimal projector orientation
    % Define the Bloch vector from the projector
    n_x = trace(P0*pauli_x);
    n_y = trace(P0*pauli_y);
    n_z = trace(P0*pauli_z);

    bloch_vector1 = [n_x,n_y,n_z];

    % Find potentially two different unitaries = set of projectors P0',
    % P1' if the system has plane or arbitrary freedom of rotation
    min_cost = cost;

    i = 1;

    % First loop to find a second Bloch vector, different from the first
    % one, with the same optimal cost
    while i<=10

        [cost2,P0_help,~] = opt_classcorr_bipartite(DensityOp,TargetQubit);

        % Define the Bloch vector from the projector
        n_x_help = trace(P0_help*pauli_x);
        n_y_help = trace(P0_help*pauli_y);
        n_z_help = trace(P0_help*pauli_z);

        bloch_vector2 = [n_x_help,n_y_help,n_z_help];

        % Calculate angle between the two vectors to see if they are the
        % same or not, if not they span a plane
        bloch_angle = real(acos((bloch_vector1*bloch_vector2')/((norm(bloch_vector1)*norm(bloch_vector2)))));

        % A better optimum has been found, restart with the new vector
        if i==1 && (min_cost-cost2)>1e-4
            bloch_vector1 = bloch_vector2;
            min_cost = min(min_cost,cost2);
            continue;
        end

        min_cost = min(min_cost,cost2);

        % Ignore optimizations which did not reach the optimal cost
        if abs(min_cost-cost2)> tolerance
            continue;
        end

        % The two optimal vectors are different, they span a plane
        if bloch_angle>0.0175 && bloch_angle < 3.140
            break;
        end

        i = i+1;
    end

    if i<=10
        % Second loop to find another potentially different Bloch vector
        j = 1;

        while j<=10

            [cost3,P0_help2,~] = opt_classcorr_bipartite(DensityOp,TargetQubit);

            % Define the Bloch vector from the projector
            n_x_help2 = trace(P0_help2*pauli_x);
            n_y_help2 = trace(P0_help2*pauli_y);
            n_z_help2 = trace(P0_help2*pauli_z);

            bloch_vector3 = [n_x_help2,n_y_help2,n_z_help2];

            % Calculate the angles between the third vector and the
            % first two vectors to see if they are the same or not, if
            % not they span more than a plane
            angle1 = real(acos((bloch_vector1*bloch_vector3')/((norm(bloch_vector1)*norm(bloch_vector3)))));
            angle2 = real(acos((bloch_vector2*bloch_vector3')/((norm(bloch_vector2)*norm(bloch_vector3)))));

            % Ignore optimizations which did not reach the optimal cost
            if abs(min_cost-cost3)> tolerance
                continue;
            end

            % The third vector differs from both other vectors
            if angle1 > 0.0175 && angle2 > 0.0175
                break;
            end

            j = j+1;
        end
    end

    % Classify the degree of freedom of the optimal orientation
    if i<=10
        % Normal vectors of the planes spanned by each pair of vectors
        normal_vector1 = cross(bloch_vector1,bloch_vector2);
        normal_vector2 = cross(bloch_vector1,bloch_vector3);
        normal_vector3 = cross(bloch_vector2,bloch_vector3);

        % Angles between the normal vectors
        norm_angle1 = real(acos((normal_vector1*normal_vector2')/((norm(normal_vector1)*norm(normal_vector2)))));
        norm_angle2 = real(acos((normal_vector1*normal_vector3')/((norm(normal_vector1)*norm(normal_vector3)))));
        norm_angle3 = real(acos((normal_vector2*normal_vector3')/((norm(normal_vector2)*norm(normal_vector3)))));

        tolerance = 1e-4;

        if (abs(norm_angle1) <= tolerance || abs(norm_angle1 - pi) <= tolerance) && ...
                (abs(norm_angle2) <= tolerance || abs(norm_angle2 - pi) <= tolerance) && ...
                (abs(norm_angle3) <= tolerance || abs(norm_angle3 - pi) <= tolerance)

            % if all the normal vectors are equivalent, then degree of
            % freedom is just one plane
            optimization_option = 'plane';

        else

            % arbitrary rotations necessary to find the correct
            % orientation
            optimization_option = 'preferred_axis';

        end

        if OptStructure == 0 && optimization_option == "plane"
            optimization_option = "none";
        end

    % No second optimal orientation has been found, the orientation is
    % unique
    elseif i>10
        optimization_option = 'none';
    end

    % Set up projectors depending on the target qubit
    if TargetQubit == 1
        P0_fullsys = kron(P0,eye(2));
        P1_fullsys = kron(P1,eye(2));
    elseif TargetQubit == 2
        P0_fullsys = kron(eye(2),P0);
        P1_fullsys = kron(eye(2),P1);
    end

    % Settings of the quasi-Newton optimizer
    Options = optimset('HessUpdate','lbfgs','TolX',1e-8,...
        'TolFun',1e-8,'MaxIter', 5000, 'display','iter', 'UseParallel',true,...
        'maxfuneval',1e10, 'LargeScale', 'off', 'GradObj', 'off', ...
        'DerivativeCheck', 'off');

    % Initialization of values for optimizations
    rot_angle_ini = pi/2;

    % Apply the projectors to the structures depending on the degree of
    % freedom of the optimal orientation
    switch optimization_option

        % Unique orientation, use the projectors as they are
        case "none"

            RhoSigmaCl = P0_fullsys*RhoSigma*P0_fullsys+P1_fullsys*RhoSigma*P1_fullsys;
            RhoMCl = P0_fullsys*RhoM*P0_fullsys+P1_fullsys*RhoM*P1_fullsys;
            RhoECl = P0_fullsys*RhoConn*P0_fullsys + P1_fullsys*RhoConn*P1_fullsys;
            RhoDCl = RhoSigma-RhoSigmaCl;

            disp("Unique orientation")

        % Plane of optimal orientations, optimize the rotation angle about
        % the normal vector of the plane
        case "plane"

            calc = @(rot_angle) cost_find_projectors(rot_angle,normal_vector1,P0_fullsys,P1_fullsys,RhoSigma,TargetQubit,2,OptStructure);
            rot_angle_opt = fminunc(calc, rot_angle_ini, Options);

            [~,RhoSigmaCl,P0_opt,P1_opt] = cost_find_projectors(rot_angle_opt,normal_vector1,P0_fullsys,P1_fullsys,RhoSigma,TargetQubit,2,OptStructure);
            RhoECl = P0_opt*RhoConn*P0_opt + P1_opt*RhoConn*P1_opt;
            RhoMCl = P0_opt*RhoM*P0_opt+P1_opt*RhoM*P1_opt;

            RhoDCl = RhoSigma-RhoSigmaCl;

            disp("Degree of freedom found: plane")

        % Arbitrary orientations, the measurement direction is chosen by
        % the user, either by the input argument or in the command window
        case "preferred_axis"

            if isempty(MeasurementDirection)
                r_axis = input("Please insert the desired measurement direction in the format [x; y; z]. Leave empty to use default direction (z-axis)");
            else
                r_axis = MeasurementDirection(:);
            end

            if numel(r_axis)< 3
                r_axis = [0;0;1];
            elseif numel(r_axis)>3
                r_axis = r_axis(1:3);
            end

            r_axis = r_axis/norm(r_axis);

            [rotated_P0, rotated_P1] = get_projectors_pref(r_axis, TargetQubit);

            RhoSigmaCl = rotated_P0*RhoSigma*rotated_P0+rotated_P1*RhoSigma*rotated_P1;
            RhoECl = rotated_P0*RhoConn*rotated_P0+rotated_P1*RhoConn*rotated_P1;
            RhoMCl = rotated_P0*RhoM*rotated_P0+rotated_P1*RhoM*rotated_P1;

            RhoDCl = RhoSigma-RhoSigmaCl;

    end

    % Reduced density operators of the two qubits after the measurement
    DeltaRho1 = partial_trace(P0_fullsys*DensityOp*P0_fullsys+P1_fullsys*DensityOp*P1_fullsys, 2);
    DeltaRho2 = partial_trace(P0_fullsys*DensityOp*P0_fullsys+P1_fullsys*DensityOp*P1_fullsys, 1);
end
