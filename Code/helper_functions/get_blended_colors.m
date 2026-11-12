function [TColors] = get_blended_colors(TBead, EBead, CBead, GridPoints, Radial)

    %!!! Compute the color data for the T-bead (total correlation) of a
    %    BEADS visualization, in which the color blends the entanglement
    %    (connected) and the compound (classical) contribution. The sign of
    %    the E-bead function values is mapped onto a yellow-blue color
    %    scale and the sign of the C-bead function values onto a red-green
    %    color scale, both with the magnitude of the T-bead function value.
    %    The two colors are then blended per vertex depending on the
    %    relative magnitude of the E- and C-bead function values: pure
    %    compound color if there is no entanglement contribution, pure
    %    entanglement color if there is no compound contribution
    %
    %    input: TBead: The function values of the T-bead (on the grid)
    %           EBead: The function values of the E-bead (on the grid)
    %           CBead: The function values of the C-bead (on the grid)
    %           GridPoints: The beads resolution [theta, phi]
    %           Radial: Shall the discrete radial-plot color scales (pure
    %                   red/green/yellow/blue with black in the center) be
    %                   used instead of the continuous color scales
    %                   (boolean)
    %
    %    output: TColors: The RGB color data of the vertices, an array of
    %                     size GridPoints(1) x GridPoints(2) x 3

    % Load the color schemes
    if Radial
        color_r = repmat([1 0 0],1000,1);
        color_g = repmat([0 1 0],1000,1);
        color_y = repmat([1 1 0],1000,1);
        color_b = repmat([0 0 1],1000,1);
        black = [0 0 0];
        color_rg = [color_g; black; color_r];
        color_yb = [color_b; black; color_y];
    else
        color_rg = flip(load('RGLinearNew20.mat').cmap);
        color_yb = flip(load('YBLinearNew20.mat').cmap);
    end

    size_yb = size(color_yb,1);

    % Get the blending angle of all vertices, 0 for pure compound and pi/2
    % for pure entanglement contribution (0/0 is treated as pi/2)
    theta_blend_matrix = atan(abs(EBead)./abs(CBead));
    theta_blend_matrix(isnan(theta_blend_matrix)) = pi/2;

    % Obtain the hue for connected from color_yb
    % range of exp values is -1 till 1
    exp_value_E = sign(EBead).*abs(TBead);

    % map exp value onto the range of 1 to size of color map
    color_index_E = round(((exp_value_E+1)/2).*size_yb+1);

    % Ensure all elements are at least 1
    color_index_E(color_index_E < 1) = 1;

    % Ensure all elements do not exceed maxIndex
    color_index_E(color_index_E > size_yb) = size_yb;

    % Look up the connected color of each vertex
    vertices = numel(color_index_E);
    Gamma_E = zeros(vertices,3);

    for i = 1:vertices
        Gamma_E(i,:) = color_yb(color_index_E(i),:);
    end

    % Obtain hue for compound
    % range of exp values is -1 till 1
    exp_value_C = sign(CBead).*abs(TBead);

    % map exp value onto the range of 1 to size of color map
    color_index_C = round(((exp_value_C+1)/2).*size_yb+1);

    % Ensure all elements are at least 1
    color_index_C(color_index_C < 1) = 1;

    % Ensure all elements do not exceed maxIndex
    color_index_C(color_index_C > size_yb) = size_yb;

    % Look up the compound color of each vertex
    vertices = numel(color_index_C);
    Gamma_C = zeros(vertices,3);

    for i = 1:vertices
        Gamma_C(i,:) = color_rg(color_index_C(i),:);
    end

    % Blend the compound and the connected color of each vertex linearly
    % with the blending angle
    Gamma_T = zeros(size(Gamma_E));

    for j = 1:vertices
        Gamma_T(j,:) = Gamma_C(j,:)+((2*theta_blend_matrix(j))/pi)*(Gamma_E(j,:)-Gamma_C(j,:));
    end

    % Reshape to the color data format of the surface plot
    TColors = reshape(Gamma_T,GridPoints(1),GridPoints(2),3);
end
