function A = matA_calvp(alpha, theta, psi, phi, omega, Q)
    %%% Function to compute matrix A for 3D vehicle sensing %%%
    % Inputs:
    %   alpha - Elevation angles (AoA)
    %   theta - Azimuth angles (AoA)
    %   psi   - Elevation angles (AoD)
    %   phi   - Azimuth angles (AoD)
    %   omega - Orientation angle
    %   Q     - Additional angle parameter for the HV orientation
    %
    % Outputs:
    %   A     - Combined matrix A

    % Number of paths (P)
    P = length(alpha);

    %%%% Matrix A computation %%%%
    % Compute components of the cosine part
    Apcos = sin(alpha) .* cos(theta) + sin(psi + Q) .* cos(phi + omega);
    Acos = [repmat(Apcos(1), P - 1, 1), diag(-Apcos(2:P))];

    % Compute components of the sine part
    Apsin = sin(alpha) .* sin(theta) + sin(psi + Q) .* sin(phi + omega);
    Asin = [repmat(Apsin(1), P - 1, 1), diag(-Apsin(2:P))];

    % Compute components of the elevation part
    Apelev = cos(alpha) + cos(psi + Q);
    Aelev = [repmat(Apelev(1), P - 1, 1), diag(-Apelev(2:P))];

    % Combine all parts to form matrix A
    A = [Acos; Asin; Aelev];
end
