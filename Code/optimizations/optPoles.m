% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [NVector, Cost, ThetaOpt, PhiOpt] = optPoles(Points)

    %!!! Optimization function for finding the pair of antipodal poles on the
    %    Bloch sphere (+n and -n) which minimizes the summed distance of a
    %    set of points to their respective nearest pole. The orientation of
    %    the pole axis is parametrized by spherical angles which are
    %    optimized starting from a random orientation
    %
    %    input: Points: The list of points on the Bloch sphere, one point
    %                   per row [x, y, z]
    %
    %    output: NVector: The optimal pole axis as Cartesian unit vector
    %            Cost: The summed distance of the points to the nearest
    %                  pole for the optimal axis
    %            ThetaOpt: The polar angle of the optimal axis, in [0, pi)
    %            PhiOpt: The azimuthal angle of the optimal axis, in
    %                    [0, 2*pi)

    % Random initial guess for the pole axis in spherical coordinates
    theta = pi*rand();
    phi = 2*pi*rand();

    optpar = [theta, phi];

    % Objective function: summed distance of the points to the nearest pole
    objFun = @(par) PoleCost(par, Points);

    % Settings of the quasi-Newton optimizer
    options = optimset('HessUpdate','lbfgs','TolX',1e-12,...
        'TolFun',1e-14,'MaxIter', 5000, 'display','iter', 'UseParallel',true,...
        'maxfuneval',1e10, 'LargeScale', 'off', 'GradObj', 'off', ...
        'DerivativeCheck', 'off');

    % Optimize the pole axis orientation
    par_opt = fminunc(objFun, optpar, options);

    % Cartesian unit vector of the optimal pole axis
    NVector = [sin(par_opt(1))*cos(par_opt(2));
               sin(par_opt(1))*sin(par_opt(2));
               cos(par_opt(1))];

    % Cost of the optimal pole axis
    Cost = PoleCost(par_opt, Points);

    % Map the optimal angles to their standard intervals
    ThetaOpt = mod(par_opt(1), pi);
    PhiOpt = mod(par_opt(2), 2*pi);
end

function Cost = PoleCost(Par, Points)

    % Cost function for a given pole axis orientation: the summed distance
    % of all points to the nearest of the two antipodal poles
    %
    %    input: Par: The spherical angles [theta, phi] of the pole axis
    %           Points: The list of points on the Bloch sphere, one point
    %                   per row
    %
    %    output: Cost: The summed distance to the nearest pole

    distances = zeros(size(Points, 1),1);

    % Cartesian unit vector of the pole axis
    n = [sin(Par(1))*cos(Par(2));
         sin(Par(1))*sin(Par(2));
         cos(Par(1))];

    % The two antipodal poles
    Q1 = n;
    Q2 = -n;

    % Distance of each point to its nearest pole
    for ii = 1:size(Points,1)
        d1 = norm(Q1-Points(ii,:)');
        d2 = norm(Q2-Points(ii,:)');
        distances(ii) = min(d1, d2);
    end

    % Minimize the total distance to the poles
    Cost = sum(distances);
end
