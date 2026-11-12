function[Cost, RhoSigmaCl, P00_, P01_, P10_, P11_] = cost_sigma_cl3(Params, RhoSigma, P00, P01, P10, P11, CorrAng, TargetQubits, OptStructure)
    
    %!!! Cost function for maximizing the overlap between tripartite M_Cl
    %    beads and the tripartite separable structure
    %
    %    input: Params: Projector initial orientations 
    %           RhoSigma: Separable structure of density operator
    %           P00/P01/P10/P11: Optimized projectors
    %           CorrAng: Are the projector orientations correlated?
    %           TargetQubits: List of target qubits
    %           OptStructure: Structure which the projectors are optimized
    %                         for (integer): 1 = separable structure (the
    %                         overlap is maximized), -1 = entangling
    %                         structure (the overlap is minimized)
    %
    %    output: cost: Overlap (multiplied by -OptStructure)
    %            RhoSigmaCl: The M_Cl component with the largest
    %                          underlying separable structure
    %            P00_/P01_/P10_,P11_: Reoptimized projectors corresponding
    %                                 to RhoSigmaCl

    % Pauli matrices
    sigma_x = [0, 1; 1, 0];
    sigma_y = [0, -1i; 1i, 0];
    sigma_z = [1, 0; 0, -1];

    % Compute rotation axis of qubit A projection in Cartesian coordinates
    nx1 = sin(Params(2)) * cos(Params(3));
    ny1 = sin(Params(2)) * sin(Params(3));
    nz1 = cos(Params(2));
    
    % Construct the corresponding matrix, n . sigma
    n_dot_sigma1 = nx1 * sigma_x + ny1 * sigma_y + nz1 * sigma_z;
    
    % Compute the local unitary rotation matrix
    U1 = cos(Params(1)/2) * eye(2) - 1i * sin(Params(1)/2) * n_dot_sigma1;

    % If the angles between the single-qubit projection directions are not
    % correlated
    if CorrAng == 0
        % Compute rotation axis of qubit B projection in Cartesian coordinates
        nx2 = sin(Params(5)) * cos(Params(6));
        ny2 = sin(Params(5)) * sin(Params(6));
        nz2 = cos(Params(5));
    
        % Construct the corresponding matrix, n . sigma
        n_dot_sigma2 = nx2 * sigma_x + ny2 * sigma_y + nz2 * sigma_z;
    
        % Compute the local unitary rotation matrix
        U2 = cos(Params(4)/2) * eye(2) - 1i * sin(Params(4)/2) * n_dot_sigma2;
    else
        U2 = U1;
    end

    % Sort the target qubit list
    targetqubits_sorted = sort(TargetQubits, 'ascend');
    % Projections of qubits 1 and 2
    if isequal(targetqubits_sorted, [1,2])
        U = kron(kron(U1, U2), eye(2));
    % Projections of qubits 1 and 3
    elseif isequal(targetqubits_sorted, [1,3])
        U = kron(kron(U1, eye(2)), U2);
    % Projections of qubits 2 and 3
    else
        U = kron(kron(eye(2), U1), U2);
    end
    
    % Compute the reoriented set of projectors
    P00_ = U*P00*U';
    P01_ = U*P01*U';
    P10_ = U*P10*U';
    P11_ = U*P11*U';

    % Extract the separable structure corresponding to the given set of
    % projectors
    RhoSigmaCl = P00_*RhoSigma*P00_ + ...
                 P01_*RhoSigma*P01_ + ...
                 P10_*RhoSigma*P10_ + ...
                 P11_*RhoSigma*P11_;
    
    % Determine the overlap with the entire classical structure, negate to
    % make it a cost function (OptStructure = 1 maximizes the overlap,
    % OptStructure = -1 minimizes it)
    Cost = -OptStructure*abs(real(trace(RhoSigmaCl * RhoSigma)));
end
    