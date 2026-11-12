function [ClassCorr, P0, P1] = opt_classcorr_bipartite(DensityOp, TargetQubit)

    %!!! Optimization function for finding the projectors of the measurement
    %    basis in which, by measuring the target qubit, the reduced density
    %    operator of the second qubit possesses the most information
    %    (minimally mixed = minimal entropy). The conditional entropy
    %    computed by cost_classcorr is minimized, i.e., the classical
    %    correlation is maximized
    %
    %    input: DensityOp: The density operator describing the state
    %           TargetQubit: The qubit which is measured
    %
    %    output: ClassCorr: The classical correlation cost value of the
    %                       optimized measurement
    %            P0/P1: The optimized single-qubit projectors on the
    %                   measured qubit

    % Depending on the dimension of the density operator, calculate the
    % number of qubits and determine the number of single qubit unitary
    % matrices which are needed
    dim = size(DensityOp,1);
    qubits = log(dim)/log(2);

    % Initialize 1 qubit unitary matrix vector for optimization
    dim_singlequbit = 2;
    U_ini = rand(1,dim_singlequbit^2)';

    % Settings of the quasi-Newton optimizer
    Options = optimset('HessUpdate','lbfgs','TolX',1e-8,...
        'TolFun',1e-8,'MaxIter', 5000, 'display','iter', 'UseParallel',true,...
        'maxfuneval',1e10, 'LargeScale', 'off', 'GradObj', 'off', ...
        'DerivativeCheck', 'off');

    % Optimize the unitary coefficients, i.e., the projector orientation
    calc = @(U_coeff) cost_classcorr(U_coeff, DensityOp, TargetQubit, dim_singlequbit, qubits);
    [U_opt, cost_final, exitflag, output] = fminunc(calc, U_ini, Options);

    % Evaluate the cost and the projectors for the optimized coefficients
    [ClassCorr,P0,P1] = cost_classcorr(U_opt, DensityOp, TargetQubit, dim_singlequbit, qubits);
end
