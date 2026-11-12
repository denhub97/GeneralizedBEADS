function [CoeffList] = traverse_and_multiply(CellArray1, CellArray2, CoeffList, DensityOp)

    %!!! Recursively traverse two nested cell arrays of identical structure,
    %    the first one containing basis operators (matrices) and the second
    %    one the corresponding scaling factors, and append the scaled
    %    expansion coefficients of an operator to a list. For each pair of
    %    leaf elements, the basis operator is multiplied with its scaling
    %    factor and the coefficient is computed as the real part of the
    %    trace of the product of the density operator with the scaled
    %    operator, i.e., the scaled expectation value
    %
    %    input: CellArray1: The (nested) cell array of the basis operators
    %           CellArray2: The (nested) cell array of the scaling factors,
    %                       must have the same structure as CellArray1
    %           CoeffList: The cell array of coefficients which is extended
    %                      (can be empty)
    %           DensityOp: The operator whose expansion coefficients are
    %                      computed
    %
    %    output: CoeffList: The cell array of coefficients, extended by the
    %                       coefficients of all leaf elements in traversal
    %                       order

    for i = 1:numel(CellArray1)
        if iscell(CellArray1{i}) && iscell(CellArray2{i})
            % If both entries are cells, recurse deeper
            [CoeffList] = traverse_and_multiply(CellArray1{i}, CellArray2{i}, CoeffList, DensityOp);
        elseif ~iscell(CellArray1{i}) && ~iscell(CellArray2{i})

            % Skip empty arrays
            if isempty(CellArray1{i}) && isempty(CellArray2{i})
                continue; % Skip this iteration
            end

            % If both entries are non-cell elements, perform multiplication
            % Scalar * Matrix
            scaled_tensor = CellArray2{i} * CellArray1{i};

            % Coefficient of the scaled operator, only consider real part
            scaled_coeff = real(trace(DensityOp*scaled_tensor));

            % Append result to the list
            CoeffList{end+1} = scaled_coeff;
        else
            % The two cell arrays must have the same structure
            error('Mismatch in structure between the two cell arrays.');
        end
    end
end
