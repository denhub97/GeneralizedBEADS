% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [NVector, Height, Cost, ThetaOpt, PhiOpt] = optTip2TipCone(Points)

    %!!! Optimization function for finding the axis orientation and height
    %    of a double cone (two cones joined tip to tip with the tips at
    %    +h*n and -h*n on the axis n through the origin) on the Bloch
    %    sphere which minimizes the summed axial distance of a set of
    %    points to their respective nearest tip. The axis orientation is
    %    parametrized by spherical angles and optimized together with the
    %    height (constrained to the interval [0, 1]), starting from a
    %    random orientation
    %
    %    input: Points: The list of points on the Bloch sphere, one point
    %                   per row [x, y, z]
    %
    %    output: NVector: The optimal cone axis as Cartesian unit vector
    %            Height: The optimal distance of the tips from the origin
    %            Cost: The summed axial distance of the points to the
    %                  nearest tip for the optimal axis and height
    %            ThetaOpt: The polar angle of the optimal axis, in [0, pi)
    %            PhiOpt: The azimuthal angle of the optimal axis, in
    %                    [0, 2*pi)

    % Random initial guess for the cone axis in spherical coordinates
    theta = pi*rand();
    phi = 2*pi*rand();

    % Initial cone height
    h = 0.5;

    optpar = [theta, phi, h];

    % Objective function: summed axial distance of the points to the
    % nearest tip
    objFun = @(par) ConeCost(par, Points);

    % Nonlinear inequality constraints: the height is restricted to the
    % interval [0, 1], no equality constraints
    nonlcon = @(par) deal([
    -par(3);        % par(3) >= 0
     par(3) - 1     % par(3) <= 1
    ], []);

    % Settings of the constrained optimizer
    options = optimset('HessUpdate','lbfgs','TolX',1e-10,...
        'TolFun',1e-10,'MaxIter', 5000, 'display','iter', 'UseParallel',true,...
        'maxfuneval',1e10, 'LargeScale', 'off', 'GradObj', 'off', ...
        'DerivativeCheck', 'off');

    % Optimize the cone axis orientation and height
    par_opt = fmincon(objFun, optpar, [], [], [], [], [], [], nonlcon, options);

    % Cartesian unit vector of the optimal cone axis
    NVector = [sin(par_opt(1))*cos(par_opt(2));
               sin(par_opt(1))*sin(par_opt(2));
               cos(par_opt(1))];

    % Optimal height of the tips
    Height = par_opt(3);

    % Cost of the optimal cone
    Cost = ConeCost(par_opt, Points);

    % Map the optimal angles to their standard intervals
    ThetaOpt = mod(par_opt(1), pi);
    PhiOpt = mod(par_opt(2), 2*pi);
end

function Cost = ConeCost(Par, Points)

    % Cost function for a given cone axis orientation and height: the
    % summed axial distance of all points to the nearest of the two tips
    %
    %    input: Par: The parameters [theta, phi, h] of the cone, i.e., the
    %                spherical angles of the axis and the tip height
    %           Points: The list of points on the Bloch sphere, one point
    %                   per row
    %
    %    output: Cost: The summed axial distance to the nearest tip

    distances = zeros(size(Points, 1),1);

    % Cartesian unit vector of the cone axis
    n = [sin(Par(1))*cos(Par(2));
         sin(Par(1))*sin(Par(2));
         cos(Par(1))];

    % The two tips on the axis
    Q1 = Par(3) * n;
    Q2 = -Par(3) * n;

    % Distance of each point to its nearest tip, measured along the axis
    for ii = 1:size(Points,1)
        d1 = abs(dot(n, Q1 - Points(ii,:)'));
        d2 = abs(dot(n, Q2 - Points(ii,:)'));
        distances(ii) = min(d1, d2);
    end

    % Minimize the total distance to the tips
    Cost = sum(distances);
end
