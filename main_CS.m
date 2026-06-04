%% Critical visibilities for 2d state families

clear
clc

% Defalut setting
opts.solver = 'mosek';
opts.verbose = 1;
opts.threads = 0;
opts.gap = 1e-6;
ell = 2; 

TotalType = 6;
rho_cell = cell(TotalType,1);
CV = zeros(TotalType,1);
run_time = zeros(TotalType,1);

for type = 1 : TotalType
    
    rho_cell{type} = generate_2d_state_family(type);
    
    t_start = tic;
    % Optimization
    
    CV(type) = Critical_visibility_canonical(rho_cell{type}, 2, opts, ell);

    run_time(type) = toc(t_start);
end

% % Without containing the information on the hierarchy level
% save('.\Num_Data\CVforD2StateFamily', "rho_cell", "CV", "run_time")

%% Critical visibilities for higher dimensional state families

clear
clc

% Defalut setting
opts.solver = 'mosek';
opts.verbose = 1;
opts.threads = 0;
opts.gap = 1e-6;
ell = 2; 

TotalType = 6;
rho_cell = cell(TotalType,1);
d  = zeros(TotalType,1);
CV = zeros(TotalType,1);
run_time = zeros(TotalType,1);

for type = 1 : TotalType % Considering six types of states.
    
    [rho_cell{type}, d(type)] = generate_quantum_sets(type); 

    t_start = tic;
    % Optimization
    
    CV(type) = Critical_visibility_aggregate(rho_cell{type}, d(type), opts, ell);
    % CV(type) = Critical_visibility_canonical(rho_cell{type}, d(type), opts, ell);

    run_time(type) = toc(t_start);

end

% % Without containing the information on the hierarchy level
% save('.\Num_Data\CVforHigherDStateFamily_aggregated', "rho_cell", "d", "CV", "run_time")

%% Classical simulation for state families

clear
clc

opts.solver = 'mosek';
opts.verbose = 1;
opts.threads = 0;
opts.gap = 1e-8;
ell = 3;

rho = generate_2d_state_family(1); d=2;
% [rho, d] = generate_quantum_sets(1);

Dtable = all_tuples(d, length(rho));
[CV, Tau_value, c_value, Pi_value, Status] = Critical_visibility_canonical(rho, d, opts, ell);

nonzeroElemIndex = find(c_value>1e-6);

save('.\Num_Data\SimulationForBB84',"rho", "d", "Dtable", ...
    "CV", "Tau_value", "c_value", "nonzeroElemIndex")

%% Witnesses for state families

clear
clc

opts.solver = 'mosek';
opts.verbose = 1;
opts.threads = 0;
opts.gap = 1e-8;
ell = 2;

rho = generate_2d_state_family(1); d=2;
% [rho, d] = generate_quantum_sets(1);
Dtable = all_tuples(d, length(rho));

[CV, W_opt] = Critical_visibility_canonical_cvx(rho, d, opts, ell);

t = 0;
for i = 1:numel(W_opt)
    rhoo = CV * rho{i} + (1-CV)/d * eye(d);
    t = t - trace(rhoo * W_opt{i});
end

save('.\Num_Data\WitnessForBB84.mat',"rho", "d", "W_opt", "t")
%% Helper function

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