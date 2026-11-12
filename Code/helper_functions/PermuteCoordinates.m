function [PermutedCoord] = PermuteCoordinates(Coordinates, PermutedPos)

    %!!! Helper function to permute coordinates for drawing coordinate
    %    axes. The elements of the coordinate cell array are rearranged
    %    according to the given permutation, e.g., to draw the axes of the
    %    coordinate system in a different orientation
    %
    %    input: Coordinates: The cell array of coordinates (e.g., one cell
    %                        entry per Cartesian axis)
    %           PermutedPos: The permutation, the ii-th entry gives the
    %                        index of the original coordinate which is
    %                        placed at position ii
    %
    %    output: PermutedCoord: The cell array of the permuted coordinates

    % Initialize the coordinate array
    PermutedCoord = cell(1,numel(PermutedPos));

    % Now rearrange the array according to the desired permutation
    for ii = 1 : numel(PermutedCoord)
        PermutedCoord{ii} = Coordinates{PermutedPos(ii)};
    end
end
