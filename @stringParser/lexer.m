function [out, varNames, pdeVarNames, eigVarNames, indVarNames] = ...
    lexer(str, problemType)
%LEXER   Lexer for string expression in CHEBFUN
%   [OUT, VARNAMES, INDVARNAME, PDEVARNAMES, EIGVARNAMES, INDVARNAMES] = 
%   LEXER(STR) 
%   Performs a lexical analysis on the string STR. 
%   Here:
%    STR:         A string of the mathematical expression we want to analyze.
%    PROBLEMTYPE: A string denoting what kind of problem we are solving. 
%                 Possible values are 'bvp', 'eig', and 'pde'.
%    OUT:         A cell array with two columns, the left is a token and the
%                 right is a label. 
%    VARNAMES:    A cell-string with names of the variables in the expression.
%    PDEVARNAMES: Contains variables that appear in PDE expressions (with a _t
%                 subscript).
%    EIGVARNAMES: The variable used to denote the eigenvalue parameter, i.e., l,
%                 lam, or lambda.
%    INDVARNAME:  A string with the name of the variable that represents the
%                 independent variable in the problem (i.e., r, t or x).
%
%   The output of this method is then passed to the LL(1) parser. For more 
%   details of compiler theory, see e.g.,
%    [1] Aho, Sethi, Ullman, Compilers: Principles, Techniques, and Tools,
%        Addison-Wesley, 1986.
%
% See also STRINGPARSER/PARSER.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Initialize an empty output.
out = [];

% A string array containing all functions which take one argument which we are
% interested in differentiating
strFun1 = char('sin', 'cos', 'tan', 'cot', 'sec', 'csc', ...
    'sinh', 'cosh', 'tanh', 'coth', 'sech', 'csch', ...
    'asin', 'acos', 'atan', 'acot', 'asec', 'acsc', ...
    'asinh', 'acosh', 'atanh', 'acoth', 'asech', 'acsch', ...
    'hypot', ...
    'asind', 'acosd', 'atand', 'acotd', 'asecd', 'acscd', 'sind', 'cosd', ...
    'tand', 'cotd', 'secd', 'cscd', ...
    'sqrt', 'exp', 'expm1', 'heaviside', 'log', 'log10', 'log2', 'log1p', ...
    'realsqrt', 'reallog', ...
    'abs', 'sign', 'var', 'std', ...
    'erf', 'erfc', 'erfcx', 'erfinv', 'erfcinv');

% String arrays containing all functions which take two or three arguments which
% we are interested in differentiating
strFun2 = char('airy', 'besselj', 'cumsum', 'diff', 'power', 'mean', ...
        'eq', 'ne', 'ge', 'gt', 'le', 'lt', 'jump');
strFun3 = char('feval', 'fred', 'volt', 'sum', 'integral');    

% Special flags to functions
strArg = char('left', 'right', 'onevar', 'start', 'end');    

% Remove all whitespace
str = strrep(str,' ','');

% Remove trailing commas
if ( strcmp(str(end), ',') )
    str(end) = [];
end

% Change quotes (") to two apostrophes ('')
str = strrep(str, '"', '''''');

% Change two minuses to one +, etc
k = 1;
while ( k < numel(str) - 1 )
    if ( strcmp(str(k), '-') )
        if ( strcmp(str(k+1), '-') )
            str(k) = '+';
            str(k+1) = [];
        elseif ( strcmp(str(k+1), '+') )
            str(k+1) = [];
        else
            k = k+1;
        end
    elseif ( strcmp(str(k), '+') )
        if ( strcmp(str(k+1), '-') )
            str(k) = '-';
            str(k+1) = [];
        elseif ( strcmp(str(k+1), '+') )
            str(k+1) = [];
        else
            k = k+1;
        end
    else
        k = k+1;
    end
end

% Vectorize the string
str = vectorize(str);

% Obtain the name of possible variables:

% List of excluded variable names
excludedNames = char(strFun1, strFun2, strFun3, strArg, 'true', 'false');

% We temporarily replace parentheses with '-' as SYMVAR will not consider
% u to be a variable in u(1).
strtmp = strrep(str, '(', '+');
strtmp = strrep(strtmp, ')', '+');
strtmp = strrep(strtmp, '''', '');
varNames = symvar(strtmp);
for k = numel(varNames):-1:1
    if ( ismember(varNames(k), excludedNames) )
        varNames(k) = [];
    end
end
varNames = unique(varNames);

% Create a cell for potential pdeVarNames and lambdaVariables
pdeVarNames = {};
eigVarNames = {};
indVarNames = {''};
pdeSubScript = [];

% r, x and t are reserved for independent variables
xLoc = strcmp('x', varNames);
varNames(xLoc) = [];
tLoc = strcmp('t', varNames);
varNames(tLoc) = [];
rLoc = strcmp('r', varNames);
varNames(rLoc) = [];

% Indicate what independent variables appear in the problems
rExists = sum(rLoc) > 0;
tExists = sum(tLoc) > 0;
xExists = sum(xLoc) > 0;
  
% Add $ to the end of the string to mark its end
str(end+1) = '$';

% In order to spot unary operators we need to store the type of the previous
% token.
prevType = 'operator';

% Loop through the input string
while ( ~strcmp(str, '$') )
    char1 = str(1);
    % Find the type of the current token
    type = myFindType(char1, prevType);
    exprEnd = 1;
    switch ( type )
        case 'num'
            % Obtain the numbers continuously (with match), their start and end
            % positions.
            regex = '[\+\-]?(([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][\+\-]?[0-9]+)?[ij]?)';
            [m, ignored, e] = regexp(str, regex, 'match', 'start', 'end');
            
            % We can run into trouble with string such as 2*3 which will become
            % 2.*3 as we vectorize. But the . here should be a part of the
            % operator, not the number, so we have to do a little check here.
            nextnum = char(m(1));
            exprEnd = e(1);
            nextChar = str(e(1)+1);
            if ( (nextnum(end) == '.') && ...
                ~isempty(regexp(nextChar,'[\+\-\*\/\^]', 'once')) )

                % Throw away the .
                nextnum(end) = [];

                % Lower the number of which we clear the string before
                exprEnd = e(1) - 1;
            end
            
            % If we encounter xi or xj, where x is a number (i.e. we have 1i,
            % 1.32i etc), we need to combine them into one token, rather than
            % lexing 1 as a number and i as a variable.
            if ( strcmp(nextChar, 'i') || strcmp(nextChar, 'j') )
                nextnum = [nextnum, nextChar];
                exprEnd = e(1) + 1; % Increase number of chars. We throw away.
            end
            out = [out ; {nextnum, 'NUM'}];
            
        case 'unary'
            % Unary operators
            exprEnd = 1;
            switch ( char1 )
                case '+'
                    out = [out ; {char1, 'UN+'}];
                case '-'
                    out = [out ; {char1, 'UN-'}];
            end
            
        case 'point'
            % If we have point, we need to check next symbol to see if we 
            % have an operator (e.g. .*) or a double (e.g. .1):
            char2 = str(2);
            type2 = myFindType(str(2), prevType);
            switch type2               
                case 'num'      % We have a floating point number
                    regex = '[0-9]+([eE][\+\-]?[0-9]+)?[ij]?';
                    [m, ignored, e] = ...
                        regexp(str, regex, 'match', 'start', 'end');

                    % Add a . and convert from cell to string
                    nextnum = ['.', char(m(1))];
                    exprEnd = e(1);
                    out = [out; {nextnum, 'NUM'}];
                    
                case 'operator' % We have a pointwise operator
                    exprEnd = 2;
                    switch char2
                        case '+'
                            out = [out ; {'.+', 'OP+'}];
                        case '-'
                            out = [out ; {'.-', 'OP-'}];
                        case '*'
                            out = [out ; {'.*', 'OP*'}];
                        case '/'
                            out = [out ; {'./', 'OP/'}];
                        case '^'
                            out = [out ; {'.^', 'OP^'}];
                    end
            end
            
        case 'operator'
            % We know that *,/ and ^ will never reach this far as we have
            % already vectorized the string. Thus, we don't have to treat
            % those operators here.
            exprEnd = 1;
            char2 = str(2);
            switch ( char1 )
                case '+'
                    out = [out ; {char1, 'OP+'}];
                case '-'
                    out = [out ; {char1, 'OP-'}];
                case '('
                    out = [out ; {char1, 'LPAR'}];
                case ')'
                    out = [out ; {char1, 'RPAR'}];
                case '='
                    if ( char2 == '=' )
                        exprEnd = 2;
                        out = [out ; {'==', 'OP=='}];
                    else
                        out = [out ; {'=', 'OP='}];                      
                    end
                case '>'
                    if ( char2 == '=' )
                        exprEnd = 2;
                        out = [out ; {'>=', 'OP>='}];
                    else
                        out = [out ; {'>', 'OP>'}];                     
                    end
                case '<'
                    if ( char2 == '=' )
                        exprEnd = 2;
                        out = [out ; {'<=', 'OP<='}];
                    else
                        out = [out ; {'<', 'OP<'}];                      
                    end
                case '~'
                    if ( char2 == '=' )
                        exprEnd = 2;
                        out = [out ; {'~=', 'OP~='}];
                    else
                        error('CHEBFUN:STRINGPARSER:lexer:UnsupportedOperator', ...
                            'Unsupported operator ~.');
                    end
            end
            
        case 'deriv'
            % The derivative symbol '
            [m, s, e] = regexp(str, '''+', 'match', 'start', 'end');
            exprEnd = e(1);
            order = e(1) - s(1) + 1;
            out = [out; {m{1}, ['DER' num2str(order)]}];
            % Find the order of the derivative
            
        case 'char'
            regex = '[a-zA-Z_][a-zA-Z_0-9]*';
            [m, ignored, e] = regexp(str, regex, 'match', 'start', 'end');
            nextString = char(m(1));   % Convert from cell to string
            exprEnd = e(1);

            % First check if we are getting pi (which should obviously be
            % treated as a number
            if ( strcmp(nextString, 'pi') )
                out = [out ; {nextString, 'NUM'}];
            % Treat l, lam and lambda specially for e-value problems
            elseif ( strcmp(problemType, 'eig') && ...
                    (strcmp(nextString, 'l') || strcmp(nextString, 'lam') || ...
                    strcmp(nextString, 'lambda')) )

                out = [out ; {nextString, 'LAMBDA'}]; %#ok<*AGROW>

                % Remove the lambda from the list of variables.
                lamMatch = cellfun(@strcmp, varNames,...
                    cellstr(repmat(nextString, length(varNames), 1)));
                varNames(lamMatch) = [];
                eigVarNames = [eigVarNames ; nextString];
            % Check if we are getting the variable we are interested in
            % differentiating w.r.t.
            elseif ( any(strcmp(nextString, varNames)) )
                % Need to treat variables which have a time derivative
                % attached on them separately.
                if ( isempty(strfind(nextString, '_')) )
                    out = [out ; {nextString, 'VAR'}];
                else
                    out = [out ; {nextString, 'PDEVAR'}];
                    % Remove varNames array, store in pdeVarNames instead.
                    pdeVarLoc = strcmp(nextString, varNames);
                    varNames(pdeVarLoc) = [];
                    pdeVarNames = [pdeVarNames ; nextString];
                    undersLoc = strfind(nextString, '_');
                    pdeSubScript = nextString(undersLoc+1:end);
                end
            % Check if this string is one of the function defined in
            % strfun1 (functions with one argument)
            elseif ( strmatch(nextString, strFun1,  'exact') )
                out = [out ; {nextString, 'FUNC1'}];     
            % Check if string is 'left' or 'right', which are args for FEVAL.
            elseif ( strmatch(nextString, strArg, 'exact') )
                out = [out; {['''' nextString ''''], 'STR'}];
            % Check if this string is one of the function defined in
            % strfun2 (functions with two arguments)
            elseif ( strmatch(nextString, strFun2, 'exact') )
                out = [out ; {nextString, 'FUNC2'}];
            % Check if this string is one of the function defined in
            % strfun2 (functions with two arguments)
            elseif ( strmatch(nextString, strFun3, 'exact') )
                out = [out ; {nextString, 'FUNC3'}];                
            % If not a function nor the variable we are interested in
            % differentiating with respect to, we treat this variable just
            % as number (this enables us e.g. to be able to differentiate w.r.t.
            % x and y separately)
            else
                out = [out ; {nextString, 'INDVAR'}];
            end
            
        case 'comma'
            out = [out ; {char1,'COMMA'}];
            
        case 'error'
            error('CHEBFUN:STRINGPARSER:lexer:unknownType', ...
                ['Invalid token ''%s'' in input.\nChebgui does not support ' ...
                '''%s'' in its input fields.'], char1, char1);
    end
    
    prevType = type;

    % Special case if we have have ) as we DON'T want unary operators then
    if ( char1 == ')' )
        prevType = 'char';
    end

    %  Throw away from string what we have already scanned
    str(1:exprEnd) = '';
end

out = [out ; {'', '$'}];

% Return the name of the independent variable. Use x if none is found.
% Check whether we have too many independent variables.
if ( strcmp(problemType, 'pde') && ((rExists + tExists + xExists) > 2) )
        error('CHEBFUN:STRINGPARSER:lexer:tooManyIndVars', ...
            'Too many independent variables in input.');
elseif ( (rExists + tExists + xExists) > 1 ) % Must be in BVP or EIG mode
        error('CHEBFUN:STRINGPARSER:lexer:tooManyIndVars', ...
            'Too many independent variables in input.');     
end

if ( rExists )
    indVarNames{1} = 'r';
elseif ( tExists )
    indVarNames{1} = 't';
elseif ( xExists )
    indVarNames{1} = 'x';
end

if ( strcmp(problemType, 'pde') )
    indVarNames{2} = pdeSubScript;
else
    indVarNames{2} = '';
end

end

function type = myFindType(str, prevtype)
% TYPE = MYFINDTYPE(STR, PREVTYPE) returns what type the current token is,
% looking at the start of the current string STR we are performing lexical
% analysis on. Here, PREVTYPE is the type of the previous token, required to
% catch unary operators.

% Change to floating point format?
% regex:  [+-]?(([0-9]+(.[0-9]*)?|.[0-9]+)([eE][+-]?[0-9]+)?)
if ( regexp(str, '[0-9]') )
    type = 'num';
% If we want to treat unary operators especially
elseif ( (strcmp(prevtype, 'operator') || strcmp(prevtype, 'unary')) && ...
        ~isempty(regexp(str,  '[+-]', 'once')) )
    type = 'unary';
elseif ( regexp(str, '[A-Za-z_]') )
    type = 'char';
elseif ( str == '.' )
    % We need to be able to distinguish between doubles and operators
    type = 'point';
elseif ( regexp(str, '\.?[\=\+\-\*\/\.\^\(\)\>\<\~]') )
    type = 'operator';
elseif ( regexp(str, '''') )
    type = 'deriv';
elseif ( strcmp(str, ',') )
    type = 'comma';
else
    % We end up here if we have encountered an unexpected type. This causes an
    % error to be thrown in the switch statement above.
    type = 'error';
end

end
