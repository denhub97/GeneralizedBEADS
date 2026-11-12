function [H] = generate_hermitian(UCoeff, RankRho)

    %!!! Generate a Hermitian matrix from a vector of real coefficients.
    %    The matrix is used as generator of a unitary, U = expm(1i*H), to
    %    parametrize unitary transformations in optimizations. The
    %    dimension of the matrix is given by the rank of the density
    %    operator (or the dimension of the single-qubit Hilbert space)
    %
    %    input: UCoeff: Vector of RankRho^2 real coefficients, the first
    %                   ones are used for the diagonal elements, the second
    %                   ones for the real and the last ones for the
    %                   imaginary parts of the upper triangular elements
    %           RankRho: The dimension of the Hermitian matrix
    %
    %    output: H: The RankRho x RankRho Hermitian matrix

    % Check that exactly RankRho^2 coefficients are provided
    assert(numel(UCoeff) == RankRho^2, ...
           "generate_hermitian:invalidCoefficients", ...
           "UCoeff must contain exactly RankRho^2 = %d coefficients, but %d were provided.", ...
           RankRho^2, numel(UCoeff));

    H = zeros(RankRho,RankRho);

    % Fill in diagonal elements
    diag_indices = 1:RankRho+1:RankRho^2;
    H(diag_indices) = UCoeff(1:RankRho);

    % Split coefficients for real and imaginary parts
    num_off_diagonal = (RankRho^2 - RankRho) / 2;
    real_coeff = UCoeff(RankRho+1:RankRho+num_off_diagonal);
    imag_coeff = UCoeff(RankRho+num_off_diagonal+1:end);

    % Fill in upper triangular elements
    upper_tri_indices = find(triu(ones(RankRho),1));
    H(upper_tri_indices) = real_coeff+1i*imag_coeff;

    % Fill in lower triangular elements by taking complex conjugate
    % diag extracts diagonal elements saves them as a column vector
    % diag(diag(H)) makes a diagonal matrix out of column vector such that
    % they can be substracted in the right way
    H = H+H'-diag(diag(H));
end
