% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [DensityOpPost, RhoPostE, RhoPostSigma] = opt_classcor_universal3(DensityOp, TargetQubits, ...
                                                                           RhoM123, RhoM12, RhoM13, RhoM23, ...
                                                                           RhoE123, RhoE12, RhoE13, RhoE23, ...
                                                                           RhoSigma123, RhoSigma12, RhoSigma13, RhoSigma23)

    %!!! Optimize the projective measurement of one or two target qubits of
    %    a tripartite state and compute the post-measurement structures of
    %    the mutual (M), the entanglement (E), and the separable (Sigma)
    %    structure of the full 3-qubit state and of its three two-qubit
    %    marginals. The projectors are optimized to minimize the
    %    conditional entropy (see opt_classcorrel3_haar_single for a single
    %    target qubit and opt_classcorrel3_haar_seq for two target qubits),
    %    taking the best result of 10 random restarts. The post-measurement
    %    structure of an operator is the sum of the operator projected with
    %    each measurement projector, sum_k P_k * Op * P_k
    %
    %    input: DensityOp: The density operator describing the tripartite
    %                      state
    %           TargetQubits: The measured qubits, either [q, q] (or a
    %                         single qubit q) for a single target qubit or
    %                         [q1, q2] for two target qubits
    %           RhoM123/RhoM12/RhoM13/RhoM23: The mutual structures of the
    %                         full state and of the marginals of qubits
    %                         1-2, 1-3, and 2-3
    %           RhoE123/RhoE12/RhoE13/RhoE23: The entanglement structures
    %                         of the full state and the marginals
    %           RhoSigma123/RhoSigma12/RhoSigma13/RhoSigma23: The separable
    %                         structures of the full state and the
    %                         marginals
    %
    %    output: DensityOpPost: Cell array with the post-measurement
    %                           structures of the mutual structures, the
    %                           entries correspond to the full state and
    %                           the marginals 1-2, 1-3, and 2-3
    %            RhoPostE: Cell array with the post-measurement
    %                      structures of the entanglement structures (same
    %                      order)
    %            RhoPostSigma: Cell array with the post-measurement
    %                          structures of the separable structures
    %                          (same order)

    % Initialize the cell arrays with the post-measurement structures
    DensityOpPost = cell(4,1);
    RhoPostE = cell(4,1);
    RhoPostSigma = cell(4,1);

    % Are we measuring a single qubit or all but a single qubit?
    if TargetQubits(1) == TargetQubits(2)
        TargetQubits = TargetQubits(1);
    end

    % Initialize the projectors and the cost of the optimization
    P0 = zeros(size(DensityOp));
    P1 = zeros(size(DensityOp));
    P00 = zeros(size(DensityOp));
    P01 = zeros(size(DensityOp));
    P10 = zeros(size(DensityOp));
    P11 = zeros(size(DensityOp));
    cost = 1e4;

    % Single target qubit
    if isscalar(TargetQubits)
        % Optimize the projectors, keep the best of 10 random restarts
        for rr = 1:10
            [cost_, P0_, P1_] = opt_classcorrel3_haar_single(DensityOp, TargetQubits);
            if cost_ < cost
                cost = cost_;
                P0 = P0_;
                P1 = P1_;
            end
        end

        % Projectors on the full system and on the two-qubit marginals,
        % the latter are traced out and normalized
        projectors = {P0, P1};
        projectors12 = {partial_trace(P0, 3)/trace(partial_trace(P0, 3)), partial_trace(P1, 3)/trace(partial_trace(P1, 3))};
        projectors13 = {partial_trace(P0, 2)/trace(partial_trace(P0, 2)), partial_trace(P1, 2)/trace(partial_trace(P1, 2))};
        projectors23 = {partial_trace(P0, 1)/trace(partial_trace(P0, 1)), partial_trace(P1, 1)/trace(partial_trace(P1, 1))};
    % Two target qubits
    else
        % Optimize the projectors, keep the best of 10 random restarts
        for rr = 1:10
            [cost_, P00_, P01_, P10_, P11_, ~, ~] = opt_classcorrel3_haar_seq(DensityOp, TargetQubits);
            if cost_ < cost
                cost = cost_;
                P00 = P00_;
                P01 = P01_;
                P10 = P10_;
                P11 = P11_;
            end
        end

        % Projectors on the full system and on the two-qubit marginals
        projectors = {P00, P01, P10, P11};
        projectors12 = {partial_trace(P00, 3), partial_trace(P01, 3), partial_trace(P10, 3), partial_trace(P11, 3)};
        projectors13 = {partial_trace(P00, 2), partial_trace(P01, 2), partial_trace(P10, 2), partial_trace(P11, 2)};
        projectors23 = {partial_trace(P00, 1), partial_trace(P01, 1), partial_trace(P10, 1), partial_trace(P11, 1)};
    end

    % Post-measurement structures of the mutual structures
    DensityOpPost{1} = post_state(RhoM123, projectors);
    DensityOpPost{2} = post_state(RhoM12, projectors12);
    DensityOpPost{3} = post_state(RhoM13, projectors13);
    DensityOpPost{4} = post_state(RhoM23, projectors23);

    % Post-measurement structures of the entanglement structures
    RhoPostE{1} = post_state(RhoE123, projectors);
    RhoPostE{2} = post_state(RhoE12, projectors12);
    RhoPostE{3} = post_state(RhoE13, projectors13);
    RhoPostE{4} = post_state(RhoE23, projectors23);

    % Post-measurement structures of the separable structures
    RhoPostSigma{1} = post_state(RhoSigma123, projectors);
    RhoPostSigma{2} = post_state(RhoSigma12, projectors12);
    RhoPostSigma{3} = post_state(RhoSigma13, projectors13);
    RhoPostSigma{4} = post_state(RhoSigma23, projectors23);

    % Discard (factor = 0) the post-measurement structure of the marginal
    % which is not used for the given set of target qubits
    factor = 0;
    if isequal(TargetQubits, 1)
        DensityOpPost{4} = factor*DensityOpPost{4};
        RhoPostE{4} = factor*RhoPostE{4};
        RhoPostSigma{4} = factor*RhoPostSigma{4};
    elseif isequal(TargetQubits, 2)
        DensityOpPost{3} = factor*DensityOpPost{3};
        RhoPostE{3} = factor*RhoPostE{3};
        RhoPostSigma{3} = factor*RhoPostSigma{3};
    elseif isequal(TargetQubits, 3)
        DensityOpPost{2} = factor*DensityOpPost{2};
        RhoPostE{2} = factor*RhoPostE{2};
        RhoPostSigma{2} = factor*RhoPostSigma{2};
    elseif isequal(TargetQubits, [1 2])
        DensityOpPost{2} = factor*DensityOpPost{2};
        RhoPostE{2} = factor*RhoPostE{2};
        RhoPostSigma{2} = factor*RhoPostSigma{2};
    elseif isequal(TargetQubits, [1 3])
        DensityOpPost{3} = factor*DensityOpPost{3};
        RhoPostE{3} = factor*RhoPostE{3};
        RhoPostSigma{3} = factor*RhoPostSigma{3};
    elseif isequal(TargetQubits, [2 3])
        DensityOpPost{4} = factor*DensityOpPost{4};
        RhoPostE{4} = factor*RhoPostE{4};
        RhoPostSigma{4} = factor*RhoPostSigma{4};
    else
        disp("Unknown set of target qubit")
    end
end

function OpPost = post_state(DensityOp, Projectors)

    %!!! Compute the post-measurement structure of an operator for a
    %    projective measurement (without resolving the outcomes), i.e.,
    %    the sum of the operator projected with each of the measurement
    %    projectors, sum_k P_k * Op * P_k
    %
    %    input: DensityOp: The density operator to be measured
    %           Projectors: Cell array of the measurement projectors P_k
    %                       (same dimension as the density operator)
    %
    %    output: OpPost: The post-measurement density operator

    % Sum the operator projected with each of the measurement projectors
    OpPost = zeros(size(DensityOp));
    for ii = 1:numel(Projectors)
        OpPost = OpPost + Projectors{ii} * DensityOp * Projectors{ii};
    end
end
