function M = fourier_bases_cross(side_dual, circ_ang)
    % Inputs:
    %   side_dual - Body-centered head direction angle (in degrees)
    %   circ_ang  - Environment/circular position angle (in degrees)
    %
    % Output:
    %   M         - Matrix of 2D Fourier cross-term basis functions

    % Convert input angles from degrees to radians
    x = side_dual / 360 * 2 * pi;
    y = circ_ang / 360 * 2 * pi;

    % Construct 2D Fourier cross-term basis up to the 2nd harmonic
    M = [cos(x).*cos(y),   sin(x).*sin(y),   cos(x).*sin(y),   sin(x).*cos(y), ...
         cos(2*x).*cos(2*y), sin(2*x).*sin(2*y), cos(2*x).*sin(2*y), sin(2*x).*cos(2*y)];
end