% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
function [DensityOp] = parse(DensityOpExpr, options)

    %!!! Parse a density operator of a 2-qubit or 3-qubit system and
    %    validate it. The density operator can be given either as a
    %    numerical matrix, which is only validated and returned, or as a
    %    text expression in Dirac notation, which is translated into a
    %    matrix and validated afterwards. In the text expression, kets
    %    |...> and bras <...| are built from the following labels, which
    %    can be concatenated to describe multiple qubits (e.g., |0+R> is
    %    the tensor product of |0>, |+>, and |R>):
    %        0, 1: computational basis
    %        +, -: Hadamard basis, (|0> +/- |1>)/sqrt(2)
    %        R, L: y-basis, (|0> +/- i|1>)/sqrt(2)
    %        Phi+, Phi-: (|00> +/- |11>)/sqrt(2)
    %        Psi+, Psi-: (|01> +/- |10>)/sqrt(2)
    %        W: (|001> + |010> + |100>)/sqrt(3)
    %        GHZ: (|000> + |111>)/sqrt(2)
    %    Besides kets and bras, the expression can contain numbers, the
    %    imaginary unit (i or j), pi, the operators + - * / ^ and
    %    parentheses, transposition ('), MATLAB matrix syntax ([1 0; 0 1]),
    %    the Kronecker product kron(ExprA, ExprB), and the functions sqrt,
    %    exp, sin, cos, tan, abs, conj, real, imag, trace, diag, eye, zeros,
    %    and ones. A combined bracket <a|b> is read as the inner product
    %    <a||b>. Other identifiers are rejected. Example:
    %    "1/2*(|0+R><0+R| + |1L-><1L-|)"
    %    If the expression evaluates to a column vector (only kets) or a
    %    row vector (only bras), it is interpreted as a pure state and cast
    %    into its projector |psi><psi|. Inside an expression, dens(ExprA)
    %    builds the projector of a column vector (|ket><ket|) or of a row
    %    vector (|ket><bra| with the ket being the conjugate transpose of
    %    the bra) from a vector expression, e.g., "0.5*dens(|Phi+>) +
    %    0.5*dens((|01>+|10>)/sqrt(2))" or "0.5*dens(|Phi+>) + 0.5*dens(<00|)". The vectors are not normalized
    %    automatically, so their norm must lead to a trace of one
    %    If the name of a parameter is specified (name-value input
    %    "Parameter"), the text expression may contain this parameter as a
    %    variable and a function handle @(Parameter) of the corresponding
    %    parametrized density operator is returned instead of a matrix, e.g.,
    %    parse("p*dens(|Phi+>) + (1-p)*eye(4)/4", "Parameter", "p"). The
    %    parameter name must be a valid MATLAB variable name and must not be
    %    R, L, or one of the identifiers listed above. As the parameter value
    %    is not known yet, the result is cast into a projector (for
    %    vectors) and validated every time the function handle is evaluated
    %
    %    input: DensityOpExpr: Either a numerical density operator (4x4 or
    %                          8x8) or a string/char array containing the
    %                          expression of the density operator
    %
    %    optional name-value input (default value in brackets):
    %           Parameter: Name of the parameter variable in the text
    %                      expression (string or char) (none)
    %
    %    output: DensityOp: The validated density operator as a matrix
    %                       (Hermitian, positive semidefinite, trace one),
    %                       or a function handle of the parameter, which
    %                       returns such a matrix, if Parameter is specified

    % Check input parameters
    arguments
        DensityOpExpr
        options.Parameter {mustBeTextScalar} = ""
    end
    Parameter = char(options.Parameter);
    assert(isnumeric(DensityOpExpr) || ischar(DensityOpExpr) || (isstring(DensityOpExpr) && isscalar(DensityOpExpr)), ...
           "parse:invalidInput", ...
           "DensityOpExpr must be a numerical matrix, a string scalar, or a char array.");

    % The parameter name must be valid
    if ~isempty(Parameter)
        check_parameter(Parameter);
    end

    % A numerical input is used as it is
    if isnumeric(DensityOpExpr)
        assert(isempty(Parameter), "parse:invalidParameter", ...
               "A parameter can only be used with a text expression, not with a numerical matrix.");
        DensityOp = DensityOpExpr;

    % A text input is translated into a matrix
    else
        Expr = char(DensityOpExpr);
        assert(isrow(Expr) && ~isempty(strtrim(Expr)), ...
               "parse:invalidInput", ...
               "The text expression must be a nonempty single line.");

        % Replace all kets and bras by references to the cell array KetCell
        [Expr, KetCell] = translate_dirac(Expr);

        % Only allow harmless characters and functions in the expression
        check_expression(Expr, Parameter);

        % A parametrized expression is returned as a function handle, which
        % is cast and validated whenever it is evaluated
        if ~isempty(Parameter)
            DensityOp = build_function_handle(Expr, KetCell, Parameter);
            return;
        end

        % Evaluate the expression
        DensityOp = evaluate_expression(Expr, KetCell);

        % A vector describes a pure state, which is cast into a density operator
        DensityOp = cast_vector(DensityOp);
    end

    % Check that the result is a valid density operator of 2 or 3 qubits
    validate_density_operator(DensityOp);
end

function [Expr, KetCell] = translate_dirac(Expr)

    %!!! Helper function to replace all kets and bras of an expression by
    %    references KetCell{k} (bras by the conjugate transpose of KetCell{k}),
    %    where KetCell contains the corresponding column vectors. Directly
    %    adjacent kets and bras (e.g., |0><1|) are joined by a multiplication
    %
    %    inputs: Expr: The text expression in Dirac notation
    %
    %    outputs: Expr: The expression with the references
    %             KetCell: Cell array of the ket vectors

    KetCell = {};
    Translated = '';
    after_ref = false;
    pos = 1;
    while pos <= numel(Expr)
        ch = Expr(pos);

        % A ket |label>
        if ch == '|'
            close_pos = find(Expr(pos+1:end) == '>', 1);
            assert(~isempty(close_pos), "parse:unclosedKet", ...
                   "The ket starting at position %d is not closed by '>'.", pos);
            KetCell{end+1} = label_to_vector(Expr(pos+1:pos+close_pos-1));
            if after_ref
                Translated = [Translated '*'];
            end
            Translated = [Translated sprintf('KetCell{%d}', numel(KetCell))];
            after_ref = true;
            pos = pos + close_pos + 1;

        % A bra <label|
        elseif ch == '<'
            close_pos = find(Expr(pos+1:end) == '|', 1);
            assert(~isempty(close_pos), "parse:unclosedBra", ...
                   "The bra starting at position %d is not closed by '|'.", pos);
            KetCell{end+1} = label_to_vector(Expr(pos+1:pos+close_pos-1));
            BraRef = [sprintf('KetCell{%d}', numel(KetCell)) ''''];
            pos = pos + close_pos + 1;
            if after_ref
                Translated = [Translated '*'];
            end
            after_ref = true;

            % A combined bracket <a|b> shares the bar between bra and ket
            Rest = Expr(pos:end);
            stop_pos = find(ismember(Rest, '<|>'), 1);
            if ~isempty(stop_pos) && stop_pos > 1 && Rest(stop_pos) == '>'
                KetCell{end+1} = label_to_vector(Rest(1:stop_pos-1));
                Translated = [Translated sprintf('(%sKetCell{%d})', [BraRef '*'], numel(KetCell))];
                pos = pos + stop_pos;
            else
                Translated = [Translated BraRef];
            end

        % A closing bracket without an opening one
        elseif ch == '>'
            error("parse:unopenedKet", "Unexpected '>' at position %d.", pos);

        % Everything else is kept
        else
            Translated = [Translated ch];
            if ~isspace(ch)
                after_ref = false;
            end
            pos = pos + 1;
        end
    end
    Expr = Translated;
end

function [Ket] = label_to_vector(Label)

    %!!! Helper function to convert the label of a ket (e.g., '0+R' or
    %    'Phi-1') into the corresponding column vector. The labels of
    %    several qubits are combined by the Kronecker product
    %
    %    inputs: Label: The text between the brackets of the ket or bra
    %
    %    outputs: Ket: The column vector of the (multi-qubit) basis state

    % Recognized labels and their vectors
    Tokens = {'0', '1', '+', '-', 'R', 'L', 'Phi+', 'Phi-', 'Psi+', 'Psi-', 'W', 'GHZ'};
    W = zeros(8,1);
    W([2 3 5]) = 1/sqrt(3);
    GHZ = zeros(8,1);
    GHZ([1 8]) = 1/sqrt(2);
    Vectors = {[1;0], [0;1], [1;1]/sqrt(2), [1;-1]/sqrt(2), [1;1i]/sqrt(2), [1;-1i]/sqrt(2), ...
               [1;0;0;1]/sqrt(2), [1;0;0;-1]/sqrt(2), [0;1;1;0]/sqrt(2), [0;1;-1;0]/sqrt(2), W, GHZ};

    % Whitespace inside the brackets is ignored
    Label = regexprep(Label, '\s', '');
    assert(~isempty(Label), "parse:emptyLabel", "A ket or bra without a label was found.");

    % Read the labels from left to right and combine them
    Ket = 1;
    pos = 1;
    while pos <= numel(Label)
        found = false;
        for t = 1:numel(Tokens)
            token_length = numel(Tokens{t});
            if strncmp(Label(pos:end), Tokens{t}, token_length)
                Ket = kron(Ket, Vectors{t});
                pos = pos + token_length;
                found = true;
                break;
            end
        end
        assert(found, "parse:unknownLabel", ...
               "Unknown basis label '%s' in the ket or bra with label '%s'. Valid labels are 0, 1, +, -, R, L, Phi+, Phi-, Psi+, Psi-, W, and GHZ.", ...
               Label(pos:end), Label);
    end
end

function [Identifiers] = allowed_identifiers()

    %!!! Helper function returning the identifiers which are allowed in a
    %    text expression (besides the parameter)
    %
    %    outputs: Identifiers: Cell array of the allowed identifiers

    Identifiers = {'KetCell', 'kron', 'sqrt', 'exp', 'sin', 'cos', 'tan', 'abs', 'conj', 'real', ...
                   'imag', 'trace', 'diag', 'eye', 'zeros', 'ones', 'pi', 'i', 'j', 'dens'};
end

function [] = check_parameter(Parameter)

    %!!! Helper function to ensure that the name of the parameter is a
    %    valid variable name which does not collide with the basis labels
    %    R and L or with the allowed identifiers
    %
    %    inputs: Parameter: The name of the parameter

    assert(isvarname(Parameter), "parse:invalidParameter", ...
           "The parameter name '%s' is not a valid variable name.", Parameter);
    assert(~ismember(Parameter, {'R', 'L'}), "parse:invalidParameter", ...
           "The parameter name must not be R or L, as these are basis labels.");
    assert(~ismember(Parameter, [allowed_identifiers(), {'Expr', 'Parameter', 'Inner'}]), "parse:invalidParameter", ...
           "The parameter name '%s' is reserved.", Parameter);
end

function [Handle] = build_function_handle(Expr, KetCell, Parameter)

    %!!! Helper function to create the function handle of a parametrized
    %    expression. Every evaluation of the handle casts vectors into
    %    projectors and validates the resulting density operator
    %
    %    inputs: Expr: The checked expression
    %            KetCell: Cell array of the ket vectors referenced in Expr
    %            Parameter: The name of the parameter in Expr
    %
    %    outputs: Handle: The function handle of the parameter

    % Function used inside the expression, captured by the handle
    dens = @dens;

    % Create the handle of the expression, it captures KetCell and dens
    try
        Inner = eval(['@(' Parameter ') ' Expr]);
    catch Exception
        error("parse:evaluationFailed", ...
              "The expression could not be converted into a function: %s", Exception.message);
    end
    Handle = @(Value) finalize_density_operator(Inner(Value));
end

function [DensityOp] = finalize_density_operator(DensityOp)

    %!!! Helper function to cast a vector into a density operator and
    %    validate the result, used for each evaluation of a parametrized
    %    expression
    %
    %    inputs: DensityOp: The evaluated expression
    %
    %    outputs: DensityOp: The validated density operator

    DensityOp = cast_vector(DensityOp);
    validate_density_operator(DensityOp);
end

function [DensityOp] = cast_vector(DensityOp)

    %!!! Helper function to cast a state vector (column or row) into its
    %    projector, other inputs are returned unchanged
    %
    %    inputs: DensityOp: The evaluated expression
    %
    %    outputs: DensityOp: The projector for a vector, otherwise the input

    if isnumeric(DensityOp) && ismatrix(DensityOp) && numel(DensityOp) > 1 && min(size(DensityOp)) == 1
        DensityOp = dens(DensityOp);
    end
end

function [] = check_expression(Expr, Parameter)

    %!!! Helper function to ensure that the expression only contains
    %    allowed characters and identifiers, as it is evaluated afterwards
    %
    %    inputs: Expr: The expression after the translation of kets and bras
    %            Parameter: The name of the parameter (empty if none)

    % Allowed characters
    assert(isempty(regexp(Expr, '[^0-9A-Za-z_\s\.\+\-\*/\^\(\)\[\]\{\},;'']', 'once')), ...
           "parse:invalidCharacter", ...
           "The expression contains an unsupported character.");

    % Allowed identifiers (numbers including exponents and imaginary suffix are removed first)
    Allowed = allowed_identifiers();
    if ~isempty(Parameter)
        Allowed{end+1} = Parameter;
    end
    Stripped = regexprep(Expr, '(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?[ij]?', '0');
    Identifiers = regexp(Stripped, '[A-Za-z_]\w*', 'match');
    Invalid = setdiff(Identifiers, Allowed);
    assert(isempty(Invalid), "parse:invalidIdentifier", ...
           "The expression contains unsupported identifiers: %s.", strjoin(Invalid, ", "));
end

function [Result] = evaluate_expression(Expr, KetCell)

    %!!! Helper function to evaluate the (checked) expression. It is a
    %    separate function such that the evaluation only sees the
    %    variables Expr and KetCell
    %
    %    inputs: Expr: The checked expression
    %            KetCell: Cell array of the ket vectors referenced in Expr
    %
    %    outputs: Result: The value of the expression

    try
        Result = eval(Expr);
    catch Exception
        error("parse:evaluationFailed", ...
              "The expression could not be evaluated: %s", Exception.message);
    end
end

function [Projector] = dens(Vector)

    %!!! Helper function to build the projector of a state vector, which
    %    can be used inside the text expressions. For a column vector (ket)
    %    it returns |ket><ket|, for a row vector (bra) the projector of the
    %    corresponding ket, i.e., |bra'><bra|
    %
    %    inputs: Vector: A column vector (ket expression) or row vector (bra
    %                    expression) with more than one element
    %
    %    outputs: Projector: The matrix of the projector

    assert(isnumeric(Vector) && ismatrix(Vector) && numel(Vector) > 1 && min(size(Vector)) == 1, ...
           "parse:invalidVector", ...
           "dens() requires a column or row vector expression (kets only or bras only).");

    % Column vector: |ket><ket|, row vector: |bra'><bra|
    if size(Vector,2) == 1
        Projector = Vector*Vector';
    else
        Projector = Vector'*Vector;
    end
end

function [] = validate_density_operator(DensityOp)

    %!!! Helper function to check that a matrix is a valid density operator
    %    of 2 or 3 qubits, i.e., a numerical, finite, square matrix of size
    %    4x4 or 8x8 which is Hermitian, positive semidefinite and has unit
    %    trace (all up to a numerical tolerance)
    %
    %    inputs: DensityOp: The matrix to be checked

    % Numerical tolerance of the checks
    tol = 1e-10;

    % Numerical, finite, two-dimensional and square
    assert(isnumeric(DensityOp) && ismatrix(DensityOp) && all(isfinite(DensityOp(:))), ...
           "parse:invalidDensityOp", ...
           "The density operator must be a finite numerical matrix.");
    assert(size(DensityOp,1) == size(DensityOp,2), ...
           "parse:invalidDimension", ...
           "The density operator must be a square matrix, but is %d x %d.", ...
           size(DensityOp,1), size(DensityOp,2));

    % Suitable number of qubits (2 or 3)
    dim = size(DensityOp,1);
    assert(dim == 4 || dim == 8, "parse:invalidDimension", ...
           "The density operator must be 4x4 or 8x8 (2 or 3 qubits), but is %d x %d. Systems with more than 3 qubits are not supported yet.", ...
           dim, dim);

    % Hermitian
    assert(norm(DensityOp - DensityOp', 'fro') <= tol*max(1, norm(DensityOp, 'fro')), ...
           "parse:notHermitian", "The density operator is not Hermitian.");

    % Unit trace
    assert(abs(trace(DensityOp) - 1) <= tol, "parse:invalidTrace", ...
           "The density operator must have trace 1, but the trace is %g.", real(trace(DensityOp)));

    % Positive semidefinite
    assert(min(eig((DensityOp + DensityOp')/2)) >= -tol, "parse:notPositive", ...
           "The density operator is not positive semidefinite.");
end
