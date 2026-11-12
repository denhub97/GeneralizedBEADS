% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [MatrixString] = Mat2QuBeads(Matrix)

    %!!! Convert a matrix into a matrix string which can be used as an
    %    input for the QuBeads software to immediately visualize the
    %    operator. The string has the format
    %    "[a11 + b11*i, a12 - b12*i; a21 + b21*i, a22 + b22*i]", i.e.,
    %    entries of a row are separated by commas and rows by semicolons,
    %    negative imaginary parts are written with a minus sign
    %
    %    input: Matrix: The (complex) matrix to be converted
    %
    %    output: MatrixString: The matrix as string in QuBeads input format

    % Get the matrix dimensions
    dim1 = size(Matrix,1);
    dim2 = size(Matrix,2);

    % Assemble the string row by row, each entry as real + imag*i
    MatrixString = "[";
    for rr = 1:dim1
        for cc = 1:dim2
            % Use a minus sign for negative imaginary parts to avoid "+ -"
            if imag(Matrix(rr,cc)) < 0
                sign_string = " - ";
            else
                sign_string = " + ";
            end
            MatrixString = append(MatrixString, ...
                           num2str(real(Matrix(rr,cc))), sign_string, ...
                           num2str(abs(imag(Matrix(rr,cc)))), "*i");

            % Separate entries within a row by comma and rows by semicolon
            if cc ~= dim2
                MatrixString = append(MatrixString,", ");
            elseif cc == dim2 && rr ~= dim1
                MatrixString = append(MatrixString,";");
            end
        end
    end
    MatrixString = append(MatrixString, "]");
end
