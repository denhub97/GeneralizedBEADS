% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [Xr, Yr, Zr] = get_filled_circle(V, r, N)

    %!!! Generate the surface coordinates of a filled circle (disk) with
    %    radius r which is perpendicular to a given vector V and centered
    %    at the origin. The circle is created in the xy plane on a polar
    %    grid and then rotated such that its normal vector points along V.
    %    The returned coordinate matrices can be plotted directly with
    %    surf()
    %
    %    input: V: The vector which is perpendicular to the circle plane
    %              (does not need to be normalized)
    %           r: The radius of the circle
    %           N: The number of grid points in radial and in angular
    %              direction
    %
    %    output: Xr/Yr/Zr: The N x N coordinate matrices of the rotated
    %                      circle

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
        % Already aligned, no rotation required
        Rmat = eye(3);
    else
        % Rotation axis and angle between the Z-axis and V
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
