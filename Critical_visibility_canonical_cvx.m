function [v_opt, W_opt, cvxStatus] = Critical_visibility_canonical_cvx(rho_cell, d, opts, ell)
%
% Functionality:
%   Compute the optimal critical visibility using a CVX-based canonical lift
%   SDP formulation for aggregated hierarchy problems. The function models and
%   solves an SDP that maximizes `v` such that the noisy states can be
%   expressed as mixtures of simulable projective outcomes.
%
% Inputs:
%   rho_cell - cell array of n density matrices (each d-by-d double)
%   d        - local dimension (positive integer)
%   opts     - (optional) struct of solver options; supported fields:
%              `solver` (default: 'mosek')
%   ell      - hierarchy level (integer >= 2)
%
% Outputs:
%   v_opt    - optimal critical visibility (scalar)
%   W_opt    - dual variables (cell array) associated with main equalities
%   cvxStatus- CVX status string returned by `cvx_status`
%
% Dependencies and Requirements (requirements):
%   - CVX (http://cvxr.com/cvx/) or YALMIP+CVX interface for modeling SDPs
%   - A solver compatible with CVX (MOSEK recommended)
%
% Example:
%   [v_opt, W_opt, status] = Critical_visibility_canonical_cvx(rho_set, 3, struct('solver','mosek'), 2);
%
% Other notes:
%   - The cvx implementation uses hermitian_semidefinite constraints and
%     constructs marginal and swap constraints for the chosen `ell`.
%   - This routine is an alternative CVX-based formulation to the
%     YALMIP implementations in the repository.
%
% Author: Mengyan-Li
% Email: mylmengyanli@gmail.com
% Date: 2026-05-06

    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'solver'), opts.solver = 'mosek'; end

    if ell < 2
        error('ell must be >= 2.');
    end

    n = numel(rho_cell);
    numMu = d^n;                 % number of deterministic strategies
    eye_d = eye(d);

    % Deterministic strategy table:
    % Dtable(mu, x) in {1,...,d}
    % Use lexicographic enumeration of tuples in {1,...,d}^n.
    Dtable = all_tuples(d, n);

    % Tuple list for hierarchy variables Pi_{i_1,...,i_ell}:
    tupleList = all_tuples(d, ell);      % size: [d^ell x ell]
    numTuples = size(tupleList, 1);
    Djoint = d^ell;

    % Precompute masks for marginal constraints:
    % idx_mask{j,s} = indices alpha such that tupleList(alpha,j) = s
    idx_mask = cell(ell, d);
    for j = 1:ell
        for s = 1:d
            idx_mask{j,s} = find(tupleList(:, j) == s);
        end
    end

    % Precompute left/right identity factors for the marginals:
    I_left  = cell(ell, 1);
    I_right = cell(ell, 1);
    for j = 1:ell
        I_left{j}  = speye(d^(j-1));
        I_right{j} = speye(d^(ell-j));
    end

    % Precompute sparse operators K{a,b} such that
    %   vec(Tr_b(V_ab * X)) = K{a,b} * vec(X)
    % with vec() in MATLAB column-major order.
    swap_ops = cell(ell, ell);
    if ell >= 2
        for a = 1:ell-1
            for b = a+1:ell
                swap_ops{a,b} = build_swap_pt_operator(d, ell, a, b);
            end
        end
    end

    cvx_begin sdp quiet
        cvx_solver(opts.solver)
        
        dual variable W{n}

        variable v
        variable c(numMu) nonnegative
        variable Tau(d, d, d, numMu) complex
        variable Pi(Djoint, Djoint, numTuples, numMu) hermitian complex
        
        maximize(v)

        subject to
            v >= 0;
            v <= 1;
            sum(c) == 1;

            % Main equalities for each x:
            % noisy POVM = deterministic mixture of Tau-blocks
            for x = 1:n
                lhs = v * rho_cell{x} + (1-v)/d * eye_d;
                rhs = cvx(zeros(d, d));
                for mu = 1:numMu
                    m = Dtable(mu, x);
                    rhs = rhs + Tau(:,:,m,mu);
                end
                W{x}:lhs == rhs ; 
            end

            % d-outcomes unit-trace projective simulability constraint
            for mu = 1:numMu
                
                % 更有助于收敛
                for m = 1:d
                    Tau(:,:,m,mu) == hermitian_semidefinite(d);
                    trace(Tau(:,:,m,mu)) == c(mu);
                end

                % sum_m Tau_{m,mu} = c_mu * I
                sum(Tau(:,:,:,mu), 3) == c(mu) * eye_d;

                % Pi PSD + trace
                for alpha = 1:numTuples
                    Pi(:,:,alpha,mu) == hermitian_semidefinite(Djoint);
                    trace(Pi(:,:,alpha,mu)) == c(mu);
                end

                % Marginal constraints:
                % sum over all indices except i_j equals
                % I ⊗ ... ⊗ Tau(i_j,mu) ⊗ ... ⊗ I
                for j = 1:ell
                    for s = 1:d
                        sumBlock = sum(Pi(:,:,idx_mask{j,s},mu), 3); 
                        rhs = kron(I_left{j}, kron(Tau(:,:,s,mu), I_right{j})); 
                        sumBlock == rhs;
                    end
                end

                % Swap constraints:
                % vec(Tr_b(V_ab * Pi_alpha,mu)) = 0 whenever tuple labels at a,b differ
                if ell >= 2
                    for a = 1:ell-1
                        for b = a+1:ell
                            for alpha = 1:numTuples
                                if tupleList(alpha, a) ~= tupleList(alpha, b)
                                    Kab = swap_ops{a,b};
                                    Kab * reshape(Pi(:,:,alpha,mu), Djoint^2, 1) == 0;
                                end
                            end
                        end
                    end
                end
            end
    cvx_end

    v_opt = v;
    % c_opt = c;
    % Tau_opt = Tau;
    % Pi_opt = Pi;
    W_opt = W;
    cvxStatus = cvx_status;
end


function tupleList = all_tuples(d, m)
% All tuples in {1,...,d}^m, lexicographic order, row-wise.
    if m == 0
        tupleList = zeros(1, 0);
        return;
    end

    grids = cell(1, m);
    [grids{:}] = ndgrid(1:d);

    tupleList = zeros(d^m, m);
    for k = 1:m
        tupleList(:, k) = grids{k}(:);
    end
end


function idx = tuple_to_index(tup, d)
% Map tup in {1,...,d}^ell to a linear index in {1,...,d^ell}
    ell = numel(tup);
    idx = 1;
    for t = 1:ell
        idx = idx + (tup(t) - 1) * d^(ell - t);
    end
end


function fullTup = insert_at(redTup, pos, val, ell)
% Insert val into position pos of a reduced tuple redTup
    fullTup = zeros(1, ell);
    p = 1;
    for t = 1:ell
        if t == pos
            fullTup(t) = val;
        else
            fullTup(t) = redTup(p);
            p = p + 1;
        end
    end
end


function K = build_swap_pt_operator(d, ell, a, b)
% Build sparse operator K such that
% vec(Tr_b(V_ab * X)) = K * vec(X),
% where X acts on (C^d)^{\otimes ell} and vec is MATLAB column-major.

    Djoint = d^ell;
    Dred   = d^(ell - 1);

    redTuples = all_tuples(d, ell - 1);

    % Each output entry has exactly d nonzero contributions
    nnz_est = (Dred^2) * d;
    rows = zeros(nnz_est, 1);
    cols = zeros(nnz_est, 1);
    vals = ones(nnz_est, 1);

    p = 0;
    for u = 1:Dred
        tupUred = redTuples(u, :);

        for v = 1:Dred
            tupVred = redTuples(v, :);
            rowIdx = u + (v - 1) * Dred;   % vec(Y) position

            for k = 1:d
                % row tuple: insert k at position b, then swap a<->b
                rowFull = insert_at(tupUred, b, k, ell);
                tmp = rowFull(a);
                rowFull(a) = rowFull(b);
                rowFull(b) = tmp;

                % column tuple: insert k at position b
                colFull = insert_at(tupVred, b, k, ell);

                r = tuple_to_index(rowFull, d);
                c = tuple_to_index(colFull, d);

                p = p + 1;
                rows(p) = rowIdx;
                cols(p) = r + (c - 1) * Djoint;   % vec(X) position
            end
        end
    end

    K = sparse(rows, cols, vals, Dred^2, Djoint^2);
end