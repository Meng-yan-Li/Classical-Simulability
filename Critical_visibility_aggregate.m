function [CV, Status] = Critical_visibility_aggregate(rho, d, opts, ell)
%
% Functionality:
%   Compute the critical visibility for a collection of input quantum states
%   using an aggregated relaxation SDP hierarchy. The function builds and solves a
%   semidefinite program to find the largest scalar `v` in [0,1] such that
%   v*rho{x} + (1-v)/d * I can be decomposed according to the classical
%   strategy constraints encoded in the SDP.
%
% Inputs:
%   rho  - cell array of n density matrices (each d-by-d double)
%   d    - local dimension (positive integer)
%   opts - (optional) struct of solver options with fields:
%          `solver`  (default: 'mosek')
%          `verbose` (default: 1)
%          `threads` (default: 0)
%          `gap`     (relative gap tolerance, default: 1e-6)
%   ell  - hierarchy level (integer >= 2)
%
% Outputs:
%   CV     - scalar in [0,1], optimal critical visibility found by the SDP
%   Status - solver information / exit status returned by YALMIP
%
% Dependencies and Requirements (requirements):
%   - YALMIP (https://yalmip.github.io/) for modeling the SDP
%   - A compatible SDP solver (MOSEK is recommended and used by default)
%
% Example:
%   rho_set = {rho1, rho2};
%   opts = struct('solver','mosek','verbose',1,'threads',4,'gap',1e-6);
%   [CV, Status] = Critical_visibility_aggregate(rho_set, 3, opts, 2);
%
% Other notes:
%   - The routine builds deterministic strategy tables and enforces marginal
%     and swap constraints for the chosen hierarchy order `ell`.
%   - Defaults are applied when `opts` is omitted or fields are missing.
%
% Author: Mengyan-Li
% Email: mylmengyanli@gmail.com
% Date: 2026-05-06

    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'solver'),  opts.solver = 'mosek'; end
    if ~isfield(opts, 'verbose'), opts.verbose = 1; end
    if ~isfield(opts, 'threads'), opts.threads = 0; end
    if ~isfield(opts, 'gap'),     opts.gap = 1e-6; end

    if ell < 2
        error('ell must be >= 2.');
    end

    n = numel(rho);
    numMu = d^n;
    eye_d = eye(d);

    % Deterministic strategy table
    Dtable = all_tuples(d, n);

    % Tuple list for hierarchy variables Pi_{i_1,...,i_ell}
    tupleList = all_tuples(d, ell);
    numTuples = size(tupleList, 1);
    Djoint = d^ell;

    % Precompute masks for marginal constraints
    idx_mask = cell(ell, d);
    for j = 1:ell
        for s = 1:d
            idx_mask{j,s} = find(tupleList(:, j) == s);
        end
    end

    % Precompute left/right identity factors
    I_left  = cell(ell, 1);
    I_right = cell(ell, 1);
    for j = 1:ell
        I_left{j}  = speye(d^(j-1));
        I_right{j} = speye(d^(ell-j));
    end

    % Precompute sparse operators K{a,b} such that
    %   vec(Tr_b(V_ab * X)) = K{a,b} * vec(X)
    swap_ops = cell(ell, ell);
    for a = 1:ell-1
        for b = a+1:ell
            swap_ops{a,b} = build_swap_pt_operator(d, ell, a, b);
        end
    end

    % ---------- YALMIP variables ----------
    Constraints = [];
    v = sdpvar(1,1);
    c = sdpvar(numMu,1);
    
    Tau_bar = cell(d, 1);
    Tau = cell(d, numMu);
    for m = 1:d    
        Tau_bar{m} = sdpvar(d, d, 'hermitian', 'complex');
        rhs = 0;
        for mu = 1:numMu
            % Structural complex Hermitian variable
            Tau{m,mu} = sdpvar(d, d, 'hermitian', 'complex');
            rhs = rhs + Tau{m,mu};
        end
        Constraints = [Constraints, Tau_bar{m}==rhs];
    end

    Pi = cell(numTuples);
    for alpha = 1:numTuples
        % Structural complex Hermitian variable
        Pi{alpha} = sdpvar(Djoint, Djoint, 'hermitian', 'complex');
    end

    mainCon = cell(n,1);

    % Basic constraints
    Constraints = [Constraints, 0 <= v <= 1];
    Constraints = [Constraints, c >= 0];
    Constraints = [Constraints, sum(c) == 1];

    % Main equalities for each x
    for x = 1:n
        lhs = v * rho{x} + (1-v)/d * eye_d;

        rhs = 0;
        for mu = 1:numMu
            m = Dtable(mu, x);
            rhs = rhs + Tau{m,mu};
        end

        mainCon{x} = (lhs == rhs);
        Constraints = [Constraints, mainCon{x}];
    end

    % Per-mu constraints
    for mu = 1:numMu

        % Tau PSD + trace + sum constraint
        sumTau = 0;
        for m = 1:d
            sumTau = sumTau + Tau{m,mu};
            Constraints = [Constraints, Tau{m,mu} >= 0];
            Constraints = [Constraints, trace(Tau{m,mu}) == c(mu)];
        end
        Constraints = [Constraints, sumTau == c(mu) * eye_d];
    end
    
    % Pi PSD + trace
    for alpha = 1:numTuples
        Constraints = [Constraints, Pi{alpha} >= 0];
        Constraints = [Constraints, trace(Pi{alpha}) == 1];
    end

    % Marginal constraints
    for j = 1:ell
        for s = 1:d
            alpha_idx = idx_mask{j,s};

            sumBlock = 0;
            for kk = 1:numel(alpha_idx)
                sumBlock = sumBlock + Pi{alpha_idx(kk)};
            end

            rhs = kron(I_left{j}, kron(Tau_bar{s}, I_right{j}));
            Constraints = [Constraints, sumBlock == rhs];
        end
    end

    % Swap constraints
    for a = 1:ell-1
        for b = a+1:ell
            Kab = swap_ops{a,b};
            for alpha = 1:numTuples
                if tupleList(alpha, a) ~= tupleList(alpha, b)
                    Constraints = [Constraints, Kab * reshape(Pi{alpha}, Djoint^2, 1) == 0];
                end
            end
        end
    end


    % ---------- Solver options ----------
    ops = sdpsettings('solver', opts.solver, 'verbose', opts.verbose, ...
                      'saveduals', 1, 'savesolverinput', 0, 'savesolveroutput', 0);

    % MOSEK options are passed through the solver-specific substructure.
    % Keep only options you need.
    ops.mosek.MSK_IPAR_NUM_THREADS = opts.threads;
    ops.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = opts.gap;

    % Optional: may reduce model size in some cases, but test carefully if you rely on duals.
    % ops.removeequalities = 1;

    % ---------- Solve ----------
    sol = optimize(Constraints, -v, ops);

    CV = value(v);
    Status = sol.info;

    % W_opt = cell(n,1);
    % for x = 1:n
    %     W_opt{x} = dual(mainCon{x});
    % end
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
                rowFull = insert_at(tupUred, b, k, ell);
                tmp = rowFull(a);
                rowFull(a) = rowFull(b);
                rowFull(b) = tmp;

                colFull = insert_at(tupVred, b, k, ell);

                r = tuple_to_index(rowFull, d);
                c = tuple_to_index(colFull, d);

                p = p + 1;
                rows(p) = rowIdx;
                cols(p) = r + (c - 1) * Djoint;
            end
        end
    end

    K = sparse(rows, cols, vals, Dred^2, Djoint^2);
end