function [Curve] = bezier(Target, ControlPoints, Radius, Color, varargin)

    %!!! bezier creates a cubic Bézier curve from 4 control points and
    %    plots it as a curved cylinder (tube) into the specified target axes
    %    object in 3D space. The curve can optionally be plotted as a dashed
    %    line
    %
    %    inputs: Target: target axes object
    %            ControlPoints: cell array of 4 control points (additional
    %                           ones will be ignored)
    %            Radius: the curved cylinder cross section radius
    %            Color: the curve color as RGB triplet
    %
    %    optional inputs (must be given together, i.e., all three or none,
    %    default: solid line):
    %            Dashed: shall the curve be plotted as a dashed line
    %                    (boolean)
    %            DashLength: the length of a dash (in number of curve
    %                        segments out of the 100 segments of the curve)
    %            DashStart: start with a dash (false) or with a blank space
    %                       (true)
    %
    %    outputs: Curve: the surface object of the Bézier curve (in case of
    %                    a dashed line, the last dash)

    if numel(varargin) > 0
        Dashed = varargin{1};
        DashLength = varargin{2};
        DashStart = varargin{3};
    else
        Dashed = false;
        DashLength = 99;
        DashStart = false;
    end

    % Hold the target axes so that the color is set correctly later
    hold(Target, 'on');
    
    % Extract control points for the cubic Bezier curve
    P0 = ControlPoints{1};
    P1 = ControlPoints{2};
    P2 = ControlPoints{3};
    P3 = ControlPoints{4};
    
    % 100 points along the curve
    t = linspace(0, 1, 100);  
    
    % Define the Bézier curve based on control points, recall that Bézier
    % curves B_n are defined:
    % B_n = sum_{i=0}^n (n choose i) (1-t)^(n-i) t^i P_i 
    % with t \in [0;1] and n = 3 for cubic Bézier curves (using n+1 control
    % points)
    B = (1 - t).^3' * P0 + ...
        3 * (1 - t).^2' .* t' * P1 + ...
        3 * (1 - t)' .* t.^2' * P2 + ...
        t'.^3 * P3;
    
    % Extract the Cartesian coordinates
    x = B(:,1);
    y = B(:,2);
    z = B(:,3);
    
    % Cylinder radius, usually dividing the input entanglement measures by
    % 25 gives a good scaling
    r_cyl = Radius/25;  

    % Resolution of circular cross-section
    n_circ = 100;

    % Precalculate a default circular cross-section in xy-plane
    theta_circ = linspace(0, 2*pi, n_circ);
    circ_x = r_cyl * cos(theta_circ);
    circ_y = r_cyl * sin(theta_circ);

    % Initialize coordinates of the curve
    [X, Y, Z] = deal([]);

    
    % Calculate the Bézier curve segmentwise
    count = 0;
    for  nn = 1:DashLength:length(x)

        if nn >= length(x)
            break;
        end
        if mod(count,2) ~= DashStart && Dashed
            count = count + 1;
            continue;
        end

        X = [];
        Y = [];
        Z = [];

        for ii = nn:nn+DashLength
            if ii >= length(x)
                break;
            end
            % Tangent vector v
            if ii < length(x)
                v = [x(ii+1) - x(ii), y(ii+1) - y(ii), z(ii+1) - z(ii)];
            else
                v = [x(ii) - x(ii-1), y(ii) - y(ii-1), z(ii) - z(ii-1)];
            end
            v = v / norm(v);
        
            % Normal vector n
            if all(abs(v - [0 0 1]) < 1e-6)
                n = [1 0 0];
            else
                n = cross(v, [0 0 1]);
                n = n / norm(n);
            end
        
            % Binormal vector b
            b = cross(v, n);
        
            % Transform the default cross section to real orientation in 3D
            % space based on tangent, normal, and binormal vector
            circ_points = x(ii)*ones(1, n_circ) + circ_x*n(1) + circ_y*b(1);
            circ_points_y = y(ii)*ones(1, n_circ) + circ_x*n(2) + circ_y*b(2);
            circ_points_z = z(ii)*ones(1, n_circ) + circ_x*n(3) + circ_y*b(3);
        
            % Add segment to coordinate containers
            X = [X; circ_points];
            Y = [Y; circ_points_y];
            Z = [Z; circ_points_z];
        end
    
        % Calculate the color data of the curve
        color_array = zeros(size(Z,1), size(Z,2), 3); 
        for rr = 1:size(Z,1)
            for cc = 1:size(Z,2)
                color_array(rr,cc,:) = Color;
            end
        end
    
        % Plot the curve
        Curve = surf(Target, X, Y, Z,'CData', color_array);
        set(Curve, 'EdgeColor', 'none');
        hold(Target, "on");
        count = count + 1;
    end
    
    
    % Still hold on for subsequent plots
    hold(Target, "on");
end
