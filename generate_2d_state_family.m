function rho_cell = generate_2d_state_family(type)
% 
% Functionality:  
%   Returns a set of 2D pure state density matrices (cell arrays) 
%   based on the input type, supporting the following predefined families:
%       1 — BB84: {|0>, |1>, |+>, |->}
%       2 — Six-state: {|0>, |1>, |+>, |->, |+i>, |-i>}
%       3 — Trine: Three equally spaced phase states
%       4 — Tetrahedral (qubit SIC): Four tetrahedral SIC pure states
% Inputs:
%   type - Scalar integer, selects the state family, values ​​1, 2, 3, or 4
% 
% Outputs:
%   rho_cell - n x 1 cell array, each element is a 2x2 density matrix (double), 
%   representing the corresponding pure state rho = |psi><psi|
% 
% Dependencies and Requirements (requirements):
%   - Only requires basic MATLAB functions (matrix operations, complex arithmetic), no additional toolboxes needed.
%   - Recommended to run on MATLAB R2016b or later for best compatibility.
%
% Example:
%   rhos = generate_2d_state_family(1); % Get the density matrix of the BB84 family
%   rho0 = rhos{1}; % First term: |0><0|
%
% Other notes:
%   The function internally normalizes the state vector to prevent numerical errors.
%
% Author: Mengyan-Li
% Email: mylmengyanli@gmail.com
% Date: 2026-05-06

    % ---------------------------
    % Computational basis
    % ---------------------------
    ket0 = [1; 0];
    ket1 = [0; 1];

    % Pauli eigenstates
    ket_plus  = (ket0 + ket1) / sqrt(2);
    ket_minus = (ket0 - ket1) / sqrt(2);
    ket_plus_i  = (ket0 + 1i * ket1) / sqrt(2);
    ket_minus_i = (ket0 - 1i * ket1) / sqrt(2);

    % ---------------------------
    % Build pure-state family
    % ---------------------------
    switch type
        case 1
            % BB84 family:
            % {|0>, |1>, |+>, |->}
            psi_cell = {ket0, ket1, ket_plus, ket_minus};

        case 2
            % Six-state family:
            % {|0>, |1>, |+>, |->, |+i>, |-i>}
            psi_cell = {ket0, ket1, ket_plus, ket_minus, ket_plus_i, ket_minus_i};

        case 3
            % Trine states:
            % |\psi_k> = (|0> + exp(2*pi*i*k/3)|1>) / sqrt(2), k = 0,1,2
            psi_cell = cell(1, 3);
            for k = 0:2
                psi_cell{k+1} = (ket0 + exp(2*pi*1i*k/3) * ket1) / sqrt(2);
            end

        case 4
            % Tetrahedral / qubit SIC states:
            % rho = (I + n_x sigma_x + n_y sigma_y + n_z sigma_z)/2
            % Here we reconstruct pure-state kets from the Bloch sphere angles.
            bloch_vecs = [ ...
                 1,  1,  1;
                 1, -1, -1;
                -1,  1, -1;
                -1, -1,  1] / sqrt(3);

            psi_cell = cell(1, 4);

            for j = 1:4
                n = bloch_vecs(j, :).';
                % Convert Bloch vector n = (sinθ cosφ, sinθ sinφ, cosθ)
                nz = max(min(real(n(3)), 1), -1);  % numerical safety
                theta_b = acos(nz);
                phi_b = atan2(real(n(2)), real(n(1)));

                psi_cell{j} = [cos(theta_b/2);
                               exp(1i*phi_b) * sin(theta_b/2)];
            end
    end
    n_states = numel(psi_cell);
    rho_cell = cell(n_states, 1);

    for j = 1:n_states
        psi = psi_cell{j};
        psi = psi / norm(psi);  % normalize for safety
        rho_cell{j} = psi * psi';
    end
end