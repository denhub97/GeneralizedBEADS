function [EoF] = get_entanglement_of_formation(DensityOp)

    %!!! Calculate the entanglement of formation of a biaprtite density operator
    %   
    %    input: DensityOp: The density operator describing the state
    %
    %    output: EoF: The entanglement of formation

    % Check input arguments
    arguments
        DensityOp (4,4) {mustBeNumeric}
    end
    
    % Compute the concurrence of the given density operator
    C = get_concurrence(DensityOp);
    
    % Compute the quantity x based on concurrence
    x = (1 + sqrt(1 - C^2)) / 2;

    % Compute the Shannon entropy function of x
    EoF = real(-x*log(x) - (1-x)*log(1-x));

    % Ensure that EoF is numerical
    if isnan(EoF)
        EoF = 0;
    end
end