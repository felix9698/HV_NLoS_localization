function B = matB_d1known(c, tdoa, psi, Q, omega, phi, d1)
    %%% Function to compute matrix B for 3D vehicle sensing %%%
    % Inputs:
    %   c     - Speed of light (or any constant factor)
    %   tdoa  - Time differences of arrival (array of P values)
    %   psi   - Elevation angles (AoD) (array of P values)
    %   Q     - Additional angle parameter for the HV orientation
    %   omega - Orientation angle
    %
    % Outputs:
    %   B     - Combined matrix B

    % Number of paths (P)
    P = length(psi);
    
    %%%% Matrix B computation %%%%
    % Compute components of matrix B
    Bcos = c * tdoa(2:P) .* sin(psi(2:P) + Q) .* cos(phi(2:P) + omega);
    Bsin = c * tdoa(2:P) .* sin(psi(2:P) + Q) .* sin(phi(2:P) + omega);
    Belev = c * tdoa(2:P) .* cos(psi(2:P) + Q);

    % compute d_1
    A1pcos = sin(psi + Q) .* cos(phi + omega) - sin(psi(1) + Q) .* cos(phi(1) + omega);
    A1psin = sin(psi + Q) .* sin(phi + omega) - sin(psi(1) + Q) .* sin(phi(1) + omega);
    A1pelev = cos(psi + Q) - cos(psi(1) + Q);

    
    % Combine all parts to form matrix B
    B = [-Bcos-A1pcos(2:P).*d1; -Bsin-A1psin(2:P).*d1; -Belev-A1pelev(2:P).*d1]; % Transpose to ensure column vectors
end
