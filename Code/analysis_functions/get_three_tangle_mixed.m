% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [ThreeTangle] = get_three_tangle_mixed(PureStates)

    %!!! Compute the convex-roof 3-tangle for a given pure state ensemble
    %    which constitutes a corresponding mixed state
    %
    %    input: PureStates: The list of pure states, the corresponding
    %                       weights need to be directly included in the
    %                       state vector, i.e., sqrt(p)*|state>
    %
    %    output: ThreeTangle: the convex-roof 3-tangle

    % Check input arguments
    arguments
        PureStates cell
    end

    % Intialize
    ThreeTangle = 0;
    
    % Compute the 3-tangle
    for ii = 1:numel(PureStates)
        % Extract the weight from the state vector
        weight = norm(PureStates{ii})^2;
        % Compute the pure state 3-tangle of the ii-th state
        [~, ~, ~, ~, tau3] = get_three_tangle((PureStates{ii}*PureStates{ii}')/weight);
        % Add to convex-roof 3-tangle
        ThreeTangle = ThreeTangle + weight * tau3;
    end
end


