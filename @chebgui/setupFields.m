function [field, allVarString, indVarName, pdeVarNames, pdeflag, ...
    eigVarNames, allVarNames]  = setupFields(guifile, input, type, allVarString)
%SETUPFIELDS   Convert input from GUI window to format useful for Chebfun.
%
% Calling sequence:
%
%   [FIELD, ALLVARSTRING, INDVARNAME, PDEVARNAMES, PDEFLAG,  EIGVARNAMES, ...
%       ALLVARNAMES] = SETUPFIELDS(GUIFILE, INPUT, TYPE, ALLVARSTRING)
%
% Here, the inputs are:
%
%   GUIFILE:        A CHEBGUI object, describing the problem shown in the GUI.
%   INPUT:          The 'String' property of a particular input field of the GUI.
%   TYPE:           The type of the field (differential equation, boundary
%                   conditions or initial guess/condition).
%   ALLVARSTRING:   A strings, containing the name of all variables that appear
%                   in a problem.
%
% The outputs are:
%   FIELD:          A string that can be converted to anonymous function,
%                   describing the INPUT variable.
%   ALLVARSTRING:   A strings, containing the name of all variables that appear
%                   in a problem.
%   INDVARNAME:     A cell-array of strings that contains the name of the
%                   independent variables that appear in the problem. The first
%                   entry corresponds to the time variable, the second to the
%                   potential time variable.
%   PDEVARNAMES:    A cell array of strings that appear on PDE format in INPUT,
%                   i.e. u_t.
%   PDEFLAG:        A vector, whose kth element is equal to 1 of the kth line of
%                   INPUT is of PDE format, e.g. u_t = u'', 0 otherwise.
%   EIGVARNAMES:    A string that indicates how the eigenvalue parameter appears
%                   in the problem, that is, either l, lam or lambda.
%   ALLVARNAMES:    A cell array of string, containing the name of all variables
%                   that appear in a problem.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% For BCs, we need to check whether varNames contains anything not found in
% varNames of the DE. Should we make the varNames of the DE as parameters? Also
% check for indVarName in deRHS.

% Number of rows in the input.
numOfRows = max(size(input));

% Binary flag for PDE detection. The kth element is going to be equal to 1 of
% the kth row of the input has PDE input in it, e.g. u_t.
pdeflag = zeros(1,numOfRows);

% What variables appear with the PDE subscript _ on them?
pdeVarNames = '';

% A temporary string for the space variable, as we don't know yet whether it is
% going to be r, t or x. This will then get replaced later.
dummys = 'DUMMYSPACE';

if ( numOfRows == 1 )
    % In this case, we only have to deal with a single line of input.
    [anFun, indVarName, allVarNames, pdeVarNames, eigVarNames, commaSeparated] = ...
        setupLine(guifile, input{1}, type);
    
    % Check whether we have too many variables involved in a PDE.
    if ( ~isempty(pdeVarNames) ) % Not a PDE (or we can't do this type yet!)
        pdeflag = true;
        
        % Are we actually working in PDE mode?
        assert(strcmpi(guifile.type, 'PDE'), ...
            'CHEBFUN:CHEBGUI:setupField:incorrectPDEmode', ...
            ['Problem specified appears to be a PDE. Please ensure you''re' ...
            'working in the correct mode in Chebgui.'])
        
        % Check whether we have a match of variable names, i.e. we want u
        % and u_t, not u and v_t:
        if ( (length(allVarNames) > 1) || (length(pdeVarNames) > 1) )
            error('CHEBFUN:CHEBGUI:setupFields:numberOfVariables', ...
                'Too many variables.');
        end
        
        % Find where the PDE variable appears:
        underScoreLocation = strfind(pdeVarNames{1},'_');
        if ( ~strcmp(allVarNames{1}, pdeVarNames{1}(1:underScoreLocation-1)) )
            error('CHEBFUN:CHEBGUI:setupFields:variableNames', ...
                'Inconsistent variable names.');
        end
    end
    
    % Only need to convert the cell to a string for DE -- information has been
    % passed for BCs
    if ( strcmp(type, 'DE') )
        allVarString = allVarNames{1};
    end
    
    % Create the string which will become the anonymous function. Put x (or t)
    % as the first argument of the anonymous function if we have a BVP or EIG.
    
    % If we are working with an eigenvalue problem, the differential equation
    % will have important information on both sides of the = sign if we are
    % solving a generalized problem. In this case, we want to keep the input in
    % two separate strings in a cell-array.
    if ( strcmp(guifile.type, 'eig') && iscell(anFun) )
        field1 = ['@(', dummys, ',', allVarString, ') ', anFun{1}];
        field2 = ['@(', dummys, ',', allVarString, ') ', anFun{2}];
        field = {field1 ; field2};
    elseif ( ~strcmp(guifile.type, 'pde') && strcmp(type, 'DE') )
        field = ['@(', dummys, ',', allVarString, ') ', anFun];
        
    % Otherwise, add variables in front of what will be anonymous functions.
    % This is not needed for 'dirichlet', 'neumann', etc... This can only happen
    % in BCs, and in that case, allVarNames will be empty
    elseif ( ~isempty(allVarNames) )
        if ( commaSeparated )
            anFun = ['[', anFun, ']'];
        end
        
        if ( strcmp(type, 'BCnew') )
            field = ['@(', dummys, ',', allVarString, ') ', anFun];
        else
            field = ['@(', allVarString, ') ', anFun, ''];
        end
    else
        field = anFun;
    end
    
else % Have a system, go through each row
    % Keep track of every variable encountered in the problem
    allVarNames = {};
    allPdeVarNames = {};
    indVarName = [];
    anFun = cell(1, numOfRows);
    % Loop through the rows
    for k = 1:numOfRows
        [anFun{k}, indVarNameTemp, varNames, pdeVarNames, eigVarNames, dummy] = ...
            setupLine(guifile, input{k}, type);
        if ( strcmp(anFun{k}([1 end]), '[]') )
            anFun{k}([1 end]) = [];
        end
        
        % Update what variables appear in the problem.
        allVarNames = [allVarNames ; varNames];
        
        % Only allow one time derivative in each line
        if ( length(pdeVarNames) > 1 )
            error('CHEBFUN:CHEBGUI:setupFields:tooManyTimeDerivatives', ...
                'Only one time derivative per line allowed')
        end
        
        if ( isempty(pdeVarNames) )
            % Indicate that no PDE variable was found in this line.
            pdeVarNames = '|';
        else
            % Indicate that this line was of PDE format.
            pdeflag(k) = 1;
        end
        
        if ( isempty(indVarName) )
            indVarName = indVarNameTemp;
        elseif ( ~myStringCompare(indVarName{1}, indVarNameTemp{1}) || ...
                ~myStringCompare(indVarName{2}, indVarNameTemp{2}) )
            error('CHEBFUN:CHEBGUI:setupFields:differentIndepVars', ...
                'Different names for independent variables detected')
        end
        
        % Update what variables appear on PDE format.
        allPdeVarNames = [allPdeVarNames ; pdeVarNames];
    end
    
    % Remove duplicate variable names
    allVarNames = unique(allVarNames);
    
    % Construct the handle part. For the DE field, we need to collect all
    % the variable names in one string. If we are working with BCs, we have
    % already passed that string in (as the parameter allVarString).
    if ( strcmp(type, 'DE') )
        allVarString = allVarNames{1};
        for varCounter = 2:length(allVarNames)
            allVarString = [allVarString, ',', allVarNames{varCounter}];
        end
    end
    
    indx = (1:numOfRows)';
    % For PDEs we need to reorder so that the order of the time derivatives
    % matches the order of the input arguments.
    if ( any(pdeflag) )
        for k = 1:numel(allPdeVarNames)
            if ( ~pdeflag(k) )
                continue % Not a PDE variable, do nothing
            end
            
            vark = allPdeVarNames{k};
            if ( strcmp(vark, '|') )
                continue  % Skip dummy pdevarname
            end
            
            % Get the kth pdevarname
            vark = vark(1:find(vark == '_', 1, 'first')-1);
            
            % Find which varname this matches
            idx3 = find(strcmp(vark,allVarNames));
            
            % Update the index list
            indx(find(indx == idx3)) = indx(k);
            indx(k) = idx3;
        end
        
        [dummy, indx] = sort(indx);     % Invert the index
        pdeflag = pdeflag(indx);        % Update the pdeflags
        allPdeVarNames(strcmp(allPdeVarNames, '|')) = []; % Delete the junk
        pdeVarNames = allPdeVarNames;
    end
    
    % Did we have PDE variables but are not actually working in PDE mode?
    assert( ~( any(pdeflag) && ~strcmpi(guifile.type, 'PDE')), ...
        'CHEBFUN:CHEBGUI:setupField:incorrectPDEmode', ...
        ['Problem specified appears to be a PDE. Please ensure that you''re' ...
        ' working in the correct mode in Chebgui.'])
    
    % If we are solving a BVP or EIG, we now need x as the first argument as
    % well. However, we don't want that variable in allVarString as we use that
    % information when setting up BCs. Create the string that goes at the start
    % of the final string.
    if ( (~strcmp(guifile.type, 'pde') && strcmp(type, 'DE')) || ...
            strcmp(type, 'BCnew') )
        fieldStart = ['@(', dummys, ',', allVarString, ') '];
    else
        fieldStart = ['@(', allVarString, ') '];
    end
    
    % Construct the function. Need to treat eig. problems separately as there we
    % can have nontrivial LHS and RHS at the same time.
    allAnFun = [];
    if ( strcmp(guifile.type, 'eig') && iscell(anFun{1}) )
        allAnFun1 = [];
        allAnFun2 = [];
        for k = 1:numOfRows
            allAnFun1 = [allAnFun1, anFun{indx(k)}{1}, ';'];
            allAnFun2 = [allAnFun2, anFun{indx(k)}{2}, ';'];
        end
        
        % Remove the last comma
        allAnFun1(end) = [];
        allAnFun2(end) = [];
        
        % Set up LHS and RHS fields
        field1 = [fieldStart, '[', allAnFun1, ']'];
        field2 = [fieldStart, '[', allAnFun2, ']'];
        field = {field1 ; field2};
    else
        % Concatenate each anonymous functions vertically.
        for k = 1:numOfRows
            allAnFun = [allAnFun, anFun{indx(k)},  '; '];
        end
        allAnFun(end-1:end) = []; % Remove the last semicolon and space.
        
        % Wrap the concatenated field with []-s.
        field = [fieldStart, '[', allAnFun, ']'];
    end
    
end

end

function [field, indVarName, varNames, pdeVarNames, eigVarNames, commaSeparated] ...
    = setupLine(guifile, input, type)
%SETUPLINE      Convert an individual input line to format useful for Chebfun.
%
% Calling sequence:
%
%   [FIELD, INDVARNAME, VARNAMES, PDEVARNAMES, EIGVARNAMES, COMMASEPARATED] =
%       SETUPLINE(GUIFILE, INPUT, TYPE)
%
% Here, the inputs are:
%
%   GUIFILE:        A CHEBGUI object, describing the problem shown in the GUI.
%   INPUT:          A single line of the 'String' property of a particular
%                   input field of the GUI.
%   TYPE:           The type of the field (differential equation, boundary
%                   conditions or initial guess/condition).
%
% The outputs are:
%   FIELD:          A string that can be converted to anonymous function,
%                   describing the INPUT variable.
%   INDVARNAME:     A cell-array of strings that contains the name of the
%                   independent variables that appear in the problem. The first
%                   entry corresponds to the time variable, the second to the
%                   potential time variable.
%   VARNAMES:       A cell array of string, containing the name of all variables
%                   that appear in the line INPUT.
%   PDEVARNAMES:    A cell array of strings that appear on PDE format in the
%                   line INPUT, i.e. u_t.
%   EIGVARNAMES:    A string that indicates how the eigenvalue parameter appears
%                   in the problem, that is, either l, lam or lambda.
%   COMMASEPARATED: Equal to 1 if the line INPUT is comma separated, for
%                   example, "u=1, v=2+x".

% Do we need to convert a BC to an anonymous function?
convertBCtoAnon = 0;

% Usually, input lines are not comma separated.
commaSeparated = 0;

if ( ~isempty(strfind(input, '@')) ) % User supplied anon. function    
    % Look at the part of the string that starts after the @(x,u,...) part.
    firstRPloc = strfind(input, ')');
    trimmedInput = input(firstRPloc+1:end);
    
    % Inspect the string, which will give us the name of the independent
    % variable (which we assume to be either r, x or t).
    [field, indVarName, varNames, pdeVarNames, eigVarNames, commaSeparated] = ...
        stringParser.str2anon(trimmedInput, guifile.type, type);
    
    return
    
elseif ( any(strcmp(type,{'BC','BCnew'})) )  % More types of syntax for BCs
    bcNum = str2num(input);
    
    % Check whether we have a number (OK), allowed strings (OK) or whether
    % we will have to convert the string to anon. function (i.e. the input
    % is on the form u' +1 = 0).
    if ( ~isempty(bcNum) )
        field = input;
        indVarName = []; % Don't need to worry about lin. func. in this case
        varNames = [];
        pdeVarNames = [];
        eigVarNames = [];
    elseif ( strcmpi(input, 'dirichlet') || strcmpi(input, 'neumann') || ...
            strcmpi(input, 'periodic') )
        % Add extra 's to allow evaluation of the string
        field = ['''', input, ''''];
        indVarName = []; % Don't need to worry about lin. func. in this case
        varNames = [];
        pdeVarNames = [];
        eigVarNames = [];
    else
        convertBCtoAnon = 1;
        guifile.type = 'bvp'; % Convert to 'BVP' type for BCs.
    end
end

if ( any(strcmp(type, {'DE', 'INIT', 'INITSCALAR'})) || convertBCtoAnon )
    % Convert to anonymous function string
    [field, indVarName, varNames, pdeVarNames, eigVarNames, commaSeparated] = ...
        stringParser.str2anon(input, guifile.type, type);
end

end

function res =  myStringCompare(str1, str2)
% Modified strcmp, if we compare with an empty string, give a match
if ( isempty(str1) || isempty(str2) )
    res = 1;
else
    res = strcmp(str1, str2);
end

end
