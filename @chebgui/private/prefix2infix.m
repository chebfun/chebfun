function [infixOut notaVAR] = prefix2infix(guifile,prefIn)

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

prefixIn = prefIn; prefCounter = 1; %#ok<NASGU> Disable warning message
NOTAVAROUT = [];
infixOut = getInfix();
notaVAR = NOTAVAROUT; NOTAVAROUT = [];

function infixOut = getInfix()
next = char(prefixIn(prefCounter,2));
if ~isempty(strmatch('OP',next))
    prefCounter = prefCounter + 1;
    exp1 = getInfix();
    exp2 = getInfix();
    % We now return different outputs depending on which operator we have.
    switch next
        case 'OP='
            infixOut = [exp1, '=', exp2];
        case 'OP+'
            if strcmp(exp1,'0') && strcmp(exp2,'0')
                infixOut = '';
            elseif  strcmp(exp1,'0')
                infixOut = exp2;
            elseif strcmp(exp2,'0')
                infixOut = exp1;
            else
                infixOut = ['(', exp1, '+', exp2, ')'];
            end
        case 'OP-'
            if strcmp(exp1,'0') && strcmp(exp2,'0')
                infixOut = '';
            elseif  strcmp(exp1,'0')
                infixOut = ['-',exp2];
            elseif strcmp(exp2,'0')
                infixOut = exp1;
            else
                infixOut = ['(', exp1, '-', exp2, ')'];
            end
        case 'OP*'
            if strcmp(exp1,'1') && strcmp(exp2,'1')
                infixOut = '1';
            elseif  strcmp(exp1,'1')
                infixOut = exp2;
            elseif strcmp(exp2,'1')
                infixOut = exp1;
            elseif  strcmp(exp1,'-1')
                infixOut = ['-',exp2];
            elseif strcmp(exp2,'-1')
                infixOut = ['-',exp1];
            elseif strcmp(exp1,'0') || strcmp(exp2,'0')
                infixOut = '0';
            elseif strcmp(exp1,'-0') || strcmp(exp2,'-0')
                infixOut = '0';
            else
                infixOut = ['(',exp1, '.*', exp2,')'];
            end
        case 'OP/'
            if strcmp(exp2,'1')
                infixOut = exp1;
            else
                infixOut = ['(',exp1, './', exp2,')'];
            end
        case 'OP^'
            infixOut = [exp1, '.^(', exp2 ,')'];
        case {'OP>','OP>=','OP<','OP<='}
            nextSym = next(3:end);
            infixOut = ['(', exp1, nextSym, exp2, ')'];
    end
elseif strcmp(next,'FUNC1')
    nextFun = char(prefixIn(prefCounter,1));
    prefCounter = prefCounter + 1;
    funcArg = getInfix();
    infixOut = [nextFun, '(', funcArg , ')'];
elseif strcmp(next,'FUNC2')
    nextFun = char(prefixIn(prefCounter,1));
    prefCounter = prefCounter + 1;
    funcArg1 = getInfix();
    funcArg2 = getInfix();
    if (strcmp(nextFun,'diff') || strcmp(nextFun,'cumsum')) && strcmp(funcArg2,'1')
        infixOut = [nextFun, '(', funcArg1, ')'];
    elseif any(strcmp(nextFun,{'fred','volt'}))
        guifiletmp = chebgui('type','bvp');
        [ignored yFredVar ignored ignored xFredVar] = lexer(guifiletmp,funcArg1);        
        anonStr = ['@(' xFredVar{1} ',' yFredVar{1} ')'];
        infixOut = [nextFun, '(', anonStr funcArg1 , ',', funcArg2 ,  ')'];
        NOTAVAROUT = [NOTAVAROUT ;  yFredVar];
    else
        infixOut = [nextFun, '(', funcArg1 , ',', funcArg2 ,  ')'];
    end
elseif strcmp(next,'FUNC3')
    nextFun = char(prefixIn(prefCounter,1));
    prefCounter = prefCounter + 1;
    funcArg1 = getInfix();
    funcArg2 = getInfix();
    funcArg3 = getInfix();
    infixOut = [nextFun, '(', funcArg1 , ',', funcArg2, ',' , funcArg3 ')'];
elseif strcmp(next(1:end-1),'DER')
    prefCounter = prefCounter + 1;
    derivArg = getInfix();
    derivOrder = next(4:end);
    if ~strcmp(derivOrder,'1')
        infixOut = ['diff(', derivArg , ',' derivOrder ')'];
    else
        infixOut = ['diff(', derivArg ')'];
    end
elseif ~isempty(strmatch('UN',next))
     nextUnary = char(prefixIn(prefCounter,1));
     prefCounter = prefCounter + 1;
     unaryArg = getInfix();
     % Only care about unary -
     if strcmp(nextUnary,'-')
         infixOut = [nextUnary, unaryArg];
     else
         infixOut = unaryArg;
     end
elseif strmatch(next,'COMMA')
    prefCounter = prefCounter + 1;
    exp1 = getInfix();
    exp2 = getInfix();
    infixOut = ['(',exp1, ',', exp2,')'];
else
    infixOut = char(prefixIn(prefCounter,1));
    prefCounter = prefCounter + 1;
end
end

newlen = length(infixOut); len = inf;
while newlen ~= len
    len = newlen;
    
    % Change two minuses to one +, etc
    k = 1;
    while k < numel(infixOut)-1
        if strcmp(infixOut(k),'-')
            if strcmp(infixOut(k+1),'-')
                infixOut(k) = '+';
                infixOut(k+1) = [];
            elseif strcmp(infixOut(k+1),'+')
                infixOut(k+1) = [];
            else
                k = k+1;
            end
        elseif strcmp(infixOut(k),'+')
            if strcmp(infixOut(k+1),'-')
                infixOut(k) = '-';
                infixOut(k+1) = [];
            elseif strcmp(infixOut(k+1),'+')
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

% We often end up with an expression on the form +(...). Prevent that from
% happening.
% if infixOut(1) == '+'
%     % Find all locations of ( and )
%     leftParLoc = strfind(infixOut,'(')
%     rightParLoc = strfind(infixOut,')')
%     if ~isempty(leftParLoc)
%         locDifference = rightParLoc - leftParLoc
%     end
% end

prefixIn = []; prefCounter = []; % Clear global variables

end