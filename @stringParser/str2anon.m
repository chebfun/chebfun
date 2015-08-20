function varargout = str2anon(str, problemType, fieldType)
%STR2ANON Converts a string on 'natural syntax form' to an anonymous function 
%         MATLAB and CHEBGUI can work with.
%
%   Calling sequence:
%       VARARGOUT = STR2ANON(STR, PROBLEMTYPE, FIELDTYPE)
%   where
%    STR:            String on 'natural syntax form'.
%    PROBLEMTYPE:    The type of problem we are solving in CHEBGUI, i.e., BVP,
%                   EIG or PDE.
%    FIELDTYPE:      What type of field of the CHEBGUI we are converting, i.e. a
%                   field for the initial guess/condition, or other fields.
%
%   If the method is called with one output argument, the output will be an
%   anonymous function on a form that is useful for Chebfun. The output will
%   start on the @(u) form of anonymous functions in Matlab
%
%   If the method is called with six output arguments, the outputs are as
%   follows:
%    anFun:          A cell array of strings, with entries corresponding to the
%                    result of parsing the input to an anonymous function form.
%                    Note that here, the strings will not start with in the @(u)
%                    form.
%    indVarNames:    A cell array of strings, with entries equal to the
%                    independent variables that appear in the input string,
%                    i.e., r, t or x.
%    varNames:       A cell arry of strings, with entries equal to the names of
%                    the dependent variables that appear in the problem, e.g., 
%                    u, v, w, ...
%    pdeVarNames:    A cell array of strings, with entries equal to the names of
%                    the variables that appear in PDE expression form, e.g. u_t.
%    eigVarNames:    A cell array of strings, with entries equal to the names of
%                    the eigenvalue parameters that appear in eigenvalue
%                    problems, i.e., 'l', 'lam' or 'lambda'.
%    commaSeparated: Equal to 1 if the input expression was comma-separated, 
%                    e.g., 'u(-1) = 0, u(1) = 1'. Equal to 0 otherwise.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    fieldType = [];
end

% Put the original string through the lexer
[lexOut, varNames, pdeVarNames, eigVarNames, indVarNames] = ...
    stringParser.lexer(str, problemType);

% Make sure we have enough variables! If parsing the initial guess, and we have
% a scalar problem, we allow that the dependent variable doesn't appear.
if ( isempty(varNames) && ~strcmp(fieldType,'INITSCALAR') )
    if ( (numel(indVarNames) > 1) && ~isempty(indVarNames{2}) )
        str = sprintf('Variables ''%s'' and ''%s'' ', indVarNames{1:2});
    else
        str = sprintf('Variable ''%s'' ', indVarNames{1});
    end

    % Throw an informative error message, depending on whether we're setting up
    % a DE/BCs or an initial guess for a system.
    if ( ~strcmp(fieldType,'INIT') )
        error('CHEBFUN:STRINGPARSER:str2anon:depvars', ...
            ['No dependent variables detected. ' str ...
             'treated as independent.']);
    else
        error('CHEBFUN:STRINGPARSER:str2anon:depvars', ...
            ['No dependent variables detected in the initial field. Input ' ...
             'must be of the form "u = x, v = 2*x, ...']);
    end
end 

% Check whether we have something like u" = u_t+v_t which we don't allow
if ( length(pdeVarNames) > 1 )
    error('CHEBFUN:STRINGPARSER:str2anon:pdeVariables', ...
        'Only one time derivative per line is allowed');
end

% Parse the output from the lexer, looking for syntax errors.
syntaxTree = stringParser.parser(lexOut);

% Convert a potential = at the top of the tree to a -. The type of the problem
% determines what complications we need to take into account, hence three
% different methods.
if ( any(strcmpi(problemType, {'bvp','ivp'})) )
    syntaxTree = stringParser.splitTree(syntaxTree);
    
    % Obtain the prefix form.
    prefixOut = stringParser.tree2prefix(syntaxTree);
    
elseif ( strcmp(problemType, 'pde') )
    [syntaxTree, pdeSign] = stringParser.splitTreePDE(syntaxTree);
    
    % Obtain the prefix form.
    prefixOut = stringParser.tree2prefix(syntaxTree);
    
    % pdeSign tells us whether we need to flip the signs. Add a unitary -
    % at the beginning of the expression
    if ( pdeSign == 1 )
        prefixOut = [{'-', 'UN-'} ; prefixOut];
    end
    
elseif ( strcmp(problemType, 'eig') )
    anFunLambda = '';
    % Convert a potential at the top of the tree = to a -.
    [syntaxTree, lambdaTree, lambdaSign] = stringParser.splitTreeEIG(syntaxTree);
    % Obtain the prefix form.
    prefixOut = stringParser.tree2prefix(syntaxTree);
    
    % If lambdaTree is not empty, convert that tree to prefix-form as well
    if ( ~isempty(lambdaTree) )
        prefixOutLambda = stringParser.tree2prefix(lambdaTree);
        
        % If the lambda part is on the LHS, we need to add a - in front of
        % the prefix expression.
        if ( lambdaSign == -1 )
            prefixOutLambda = [{'-', 'UN-'} ; prefixOutLambda];
        end
        
        % If we're in EIG mode, we want to replace lambda by 1
        if ( ~isempty(eigVarNames) )
            eigvarLoc = find(ismember(prefixOutLambda(:,2), 'LAMBDA') == 1);
            prefixOutLambda(eigvarLoc,1) = ...
                cellstr(repmat('1', length(eigvarLoc), 1));
            prefixOutLambda(eigvarLoc,2) = ...
                cellstr(repmat('NUM', length(eigvarLoc), 1));
        end

        % Change it to infix form and remove uneccessary parenthesis.
        infixOutLambda = stringParser.pref2inf(prefixOutLambda);
        anFunLambda = stringParser.parSimp(infixOutLambda);
    end
    
end

% Check whether we have equations divided by commas. This will only have
% happened if we have any commas left in the prefix expression
commaSeparated = any(strcmp(prefixOut(:,2), 'COMMA'));

% Convert the prefix form to infix form.
[infixOut, notaVAR] = stringParser.pref2inf(prefixOut);

% Get rid of unnecessary parenthesis.
anFun = stringParser.parSimp(infixOut);

% Remove misinterpreted VARs (from Fred, Volt, etc)
for k = numel(varNames):-1:1
    if ( any(strcmp(varNames{k}, notaVAR)) )
        varNames(k) = [];
    end
end 
 
% Convert the cell array varNames into one string. Not required when we're
% working with the initial guess of scalar problems. If varNames is empty
% for other kind of problems, an error would already have been thrown.
if ( ~isempty(varNames) && ~strcmp(fieldType, 'INITSCALAR') )
    varString = varNames{1};
    for varCounter = 2:length(varNames)
        varString = [varString, ',', varNames{varCounter}]; %#ok<AGROW>
    end

    if ( length(varNames) == 1 )
        anFunComplete = ['@(' varString ') ' anFun];
    else
        anFunComplete = ['@(' varString ') [' anFun ']'];
    end
end

% Also return the lambda part if we are in EIG mode
if ( strcmp(problemType, 'eig') && ~isempty(anFunLambda) )
    anFunLambdaComplete = ['@(' varString ') ' anFunLambda];
    anFunComplete = {anFunComplete ; anFunLambdaComplete};
    anFun = {anFun ; anFunLambda};
end

% The output depends on the number of output variables.
if ( nargout == 1 )
    varargout{1} = anFunComplete;
else
    varargout = {anFun, indVarNames, varNames, ...
                 pdeVarNames, eigVarNames, commaSeparated};
end
                 

end
