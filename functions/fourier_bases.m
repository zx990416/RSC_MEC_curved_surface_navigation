function M = fourier_bases(side_dual)
    % Input:
    %   side_dual - Angle vector in degrees
    %
    % Output:
    %   M         - Matrix of 1D Fourier basis functions up to 4th harmonic

    % Convert input angle from degrees to radians
    side_dual = side_dual / 360 * 2 * pi;

    % Construct Fourier basis matrix
    M = [cos(side_dual),   sin(side_dual), ...
         cos(2*side_dual), sin(2*side_dual), ...
         cos(3*side_dual), sin(3*side_dual), ...
         cos(4*side_dual), sin(4*side_dual)];
end