function [rho_cell, d] = generate_quantum_sets(type)
%
% Functionality:
%   Generate several predefined sets of quantum states (density matrices)
%   for testing and demonstrations. Supported `type` values select between
%   qutrit MUB sets, SIC states, hybrid classical mixtures, and random
%   ensembles.
%
% Inputs:
%   type - scalar integer selecting the state set. Supported values:
%          1,2,3 : sets of 2,3,4 MUBs in d=3 (qutrit MUBs)
%          4     : 3-dimensional SIC pure states (9 states)
%          5,6   : hybrid diagonal states in d=4 or d=5
%          7     : random density matrices in d=5 (requires numStates variable)
%
% Outputs:
%   rho_cell - cell array of density matrices (each d-by-d double)
%   d        - local Hilbert space dimension
%
% Dependencies and Requirements (requirements):
%   - Only MATLAB base functions are required. No extra toolboxes needed.
%
% Example:
%   [rhos, d] = generate_quantum_sets(1); % get qutrit MUB states
%
% Other notes:
%   - For `type==7` the variable `numStates` must be defined in the caller
%     workspace or the function should be adapted to accept an explicit
%     number of random states.
%
% Author: Mengyan-Li
% Email: mylmengyanli@gmail.com
% Date: 2026-05-06


    switch type

        case {1, 2, 3}
            % 3-dimensional mutually unbiased bases, MUBs
            d = 3;
            numBases = type + 1;   % type=1 -> 2 bases, type=2 -> 3 bases, type=3 -> 4 bases
            rho_cell = {};
            mub_cell = generate_qutrit_mubs();
            idx = 1;
            for b = 1:numBases
                B = mub_cell{b};
                for k = 1:d
                    psi = B(:, k);
                    rho_cell{idx, 1} = psi * psi';
                    idx = idx + 1;
                end
            end

        case 4
            % 3-dimensional SIC-POVM pure states, total 9
            d = 3;
            omega = exp(2*pi*1i/3);
            V = [
                0,      0,        0,       -1,     -omega, -omega^2,  1,       1,        1;
                1,      1,        1,        0,      0,      0,       -1,      -omega,   -omega^2;
               -1,     -omega,   -omega^2, 1,      1,      1,        0,       0,        0
            ] / sqrt(2);
            rho_cell = cell(9, 1);
            for k = 1:9
                psi = V(:, k);
                rho_cell{k} = psi * psi';
            end

        case {5,6}
            % 4-dimensional hybrid state 
            % {|0><0|,|1><1|,|2><2|,(|0><0|+|1><1|+|2><2|+|3><3|)/4}
            if type == 5
                d = 4;
            else
                d = 5;
            end
            basis = eye(d);
            rho_cell = cell(d, 1);
            for k = 1:d-1
                psi = basis(:, k);
                rho_cell{k} = psi * psi';
            end
            rho_cell{d} = ones(d,1) * ones(d,1)' / d;

        case 7
            % Random density matrices in dimension 5
            d = 5;
            rho_cell = cell(numStates, 1);
            for k = 1:numStates
                A = randn(d) + 1i * randn(d);
                rho = A * A';
                rho = rho / trace(rho);
                rho_cell{k} = rho;
            end

        otherwise
            error('Unsupported type. type must be one of 1, 2, 3, 4, 5, 6, 7.');
    end

end


function mub_cell = generate_qutrit_mubs()
%GENERATE_QUTRIT_MUBS Generate a complete set of 4 MUBs in dimension 3.

    d = 3;
    omega = exp(2*pi*1i/d);
    n = (0:d-1).' ;

    mub_cell = cell(4, 1);

    % Basis 1: computational basis
    mub_cell{1} = eye(d);

    % Bases 2-4: quadratic construction for prime dimension d=3
    for b = 0:d-1
        B = zeros(d, d);
        for m = 0:d-1
            psi = omega .^ (b * n.^2 + m * n);
            B(:, m+1) = psi / sqrt(d);
        end
        mub_cell{b+2} = B;
    end

end
