function [infixOut, notaVAR] = pref2inf(prefixIn)
%PREF2INF   Convert an expression on prefix form to infix form
%   [INFIXOUT, NOTAVAR] = PREF2INF(PREFIXIN) goes recursively through the
%   expression PREFIXIN, which is on prefix form. The output, INFIXOUT, is a
%   string, representing the expression on infix form. NOTAVAR is a string, that
%   corresponds to the second variable in the kernel for FRED and VOLT
%   operators.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

prefCounter = 1;
NOTAVAROUT = [];
infixOut = getInfix();
notaVAR = NOTAVAROUT;
NOTAVAROUT = [];

function infixOut = getInfix()
%GETINFIX   Return the infix form.

next = char(prefixIn(prefCounter,2));

% Steep through the prefix expression, dealing with each different type of
% tokens separately.
if ( ~isempty(strmatch('OP', next)) )
    % Increase the counter.
    prefCounter = prefCounter + 1;
    
    % We are dealing with a binary operator, so obtain two infix expressions.
    exp1 = getInfix();
    exp2 = getInfix();
    
    % We now return different outputs depending on which operator we have.
    switch next
        case 'OP='
            infixOut = [exp1, '=', exp2];
            
        case 'OP+'
            if ( strcmp(exp1, '0') && strcmp(exp2, '0') )
                % Don't want to introduce unnecessary zeros.
                infixOut = '';
            elseif  ( strcmp(exp1, '0') )
                infixOut = exp2;
            elseif ( strcmp(exp2, '0') )
                infixOut = exp1;
            else
                infixOut = ['(', exp1, '+', exp2, ')'];
            end
            
        case 'OP-'
            if ( strcmp(exp1, '0') && strcmp(exp2, '0') )
                infixOut = '';
            elseif  ( strcmp(exp1, '0') )
                infixOut = ['-',exp2];
            elseif ( strcmp(exp2, '0') )
                infixOut = exp1;
            else
                infixOut = ['(', exp1, '-', exp2, ')'];
            end
            
        case 'OP*'
            % Check whether we have some special arguments, which we know we can
            % simplify.
            if ( strcmp(exp1, '1') && strcmp(exp2, '1') )
                infixOut = '1';
            elseif  ( strcmp(exp1, '1') )
                infixOut = exp2;
            elseif ( strcmp(exp2, '1') )
                infixOut = exp1;
            elseif  ( strcmp(exp1, '-1') )
                infixOut = ['-', exp2];
            elseif ( strcmp(exp2, '-1') )
                infixOut = ['-', exp1];
            elseif ( strcmp(exp1, '0') || strcmp(exp2, '0') )
                infixOut = '0';
            elseif ( strcmp(exp1, '-0') || strcmp(exp2, '-0') )
                infixOut = '0';
            else
                infixOut = ['(', exp1, '.*', exp2, ')'];
            end
            
        case 'OP/'
            if ( strcmp(exp2, '1') )
                % Division by 1 is trivial.
                infixOut = exp1;
            else
                infixOut = ['(', exp1, './', exp2, ')'];
            end
            
        case 'OP^'
            infixOut = [exp1, '.^(', exp2, ')'];
            
        case {'OP>', 'OP>=', 'OP<', 'OP<='}
            nextSym = next(3:end);
            infixOut = ['(', exp1, nextSym, exp2, ')'];
    end
    
elseif ( strcmp(next, 'FUNC1') )
    % A function which has one argument.
    nextFun = char(prefixIn(prefCounter, 1));
    prefCounter = prefCounter + 1;
    funcArg = getInfix();
    infixOut = [nextFun, '(', funcArg , ')'];
    
elseif ( strcmp(next, 'FUNC2') )
    % A function with two arguments.
    nextFun = char(prefixIn(prefCounter, 1));
    prefCounter = prefCounter + 1;
    % Obtain two infix expressions.
    funcArg1 = getInfix();
    funcArg2 = getInfix();
    
    % Some methods need a special treatment.
    if ( (strcmp(nextFun, 'diff') || strcmp(nextFun, 'cumsum')) && ...
            strcmp(funcArg2, '1') )
        infixOut = [nextFun, '(', funcArg1, ')'];
        
    elseif ( any(strcmp(nextFun, {'fred', 'volt'})) )
        % For fred() and volt(), we need to find what arguments denote the
        % kernel variable.
        [ignored1, yFredVar, ignored2, ignored3, xFredVar] = ...
            stringParser.lexer(funcArg1, 'bvp');
        anonStr = ['@(' xFredVar{1} ',' yFredVar{1} ')'];
        infixOut = [nextFun, '(', anonStr, funcArg1 , ',', funcArg2 ,  ')'];
        NOTAVAROUT = [NOTAVAROUT ;  yFredVar];
        
    else
        infixOut = [nextFun, '(', funcArg1 , ',', funcArg2 ,  ')'];
    end
    
elseif ( strcmp(next, 'FUNC3') )
    % Methods with three arguments.
    nextFun = char(prefixIn(prefCounter,1));
    prefCounter = prefCounter + 1;
    funcArg1 = getInfix();
    funcArg2 = getInfix();
    funcArg3 = getInfix();
    infixOut = [nextFun, '(', funcArg1 , ',', funcArg2, ',' , funcArg3 ')'];
    
elseif ( strcmp(next(1:end-1), 'DER') )
    % Differentiation via the ' symbol.
    prefCounter = prefCounter + 1;
    derivArg = getInfix();
    derivOrder = next(4:end);
    if ( ~strcmp(derivOrder, '1') )
        infixOut = ['diff(', derivArg, ',' derivOrder, ')'];
    else
        infixOut = ['diff(', derivArg, ')'];
    end
    
elseif ( ~isempty(strmatch('UN', next)) )
    % Unary operators
    nextUnary = char(prefixIn(prefCounter, 1));
    prefCounter = prefCounter + 1;
    unaryArg = getInfix();
    % Only care about unary -
    if ( strcmp(nextUnary, '-') )
        infixOut = [nextUnary, unaryArg];
    else
        infixOut = unaryArg;
    end
    
elseif ( strmatch(next, 'COMMA') )
    % The , symbol.
    prefCounter = prefCounter + 1;
    exp1 = getInfix();
    exp2 = getInfix();
    infixOut = ['(', exp1, ';', exp2, ')'];
else
    % We have arrived at a terminal token. Return it and increase the counter.
    infixOut = char(prefixIn(prefCounter,1));
    prefCounter = prefCounter + 1;
end

end

% Initialization for simplification.
newlen = length(infixOut);
len = inf;

% Do some simplifications on the returned string. In particular, we want to
% change -- to +, a -+ to -, a +- to -, and ++ to +.
while ( newlen ~= len )     % Loop until length doesn't change anymore.
    len = newlen;
    
    k = 1;
    while ( k < (numel(infixOut) - 1) )
        if ( strcmp(infixOut(k), '-') )
            if ( strcmp(infixOut(k+1), '-') )
                infixOut(k) = '+';
                infixOut(k+1) = [];
            elseif ( strcmp(infixOut(k+1), '+') )
                infixOut(k+1) = [];
            else
                k = k+1;
            end
        elseif ( strcmp(infixOut(k), '+') )
            if ( strcmp(infixOut(k+1), '-') )
                infixOut(k) = '-';
                infixOut(k+1) = [];
            elseif ( strcmp(infixOut(k+1), '+') )
                infixOut(k+1) = [];
            else
                k = k+1;
            end
        else
            k = k+1;
        end
    end
        
    newlen = length(infixOut);

end

% Clear global variables
prefixIn = [];
prefCounter = [];

end
