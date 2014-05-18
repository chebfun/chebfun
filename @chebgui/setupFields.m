function [field, allVarString, indVarName, pdeVarNames, pdeflag, eigVarNames, ...
    allVarNames]  = setupFields(guifile, input, type, allVarString)

% TODO:  Documentation.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/chebfun/ for Chebfun information.


% For BCs, we need to check whether varNames contains anything not found in
% varNames of the DE. Should we make the varNames of the DE as parameters?
% Setja DE varNames sem parametra? Also check for indVarName in deRHS.

numOfRows = max(size(input)); %numel(input,1);
pdeflag = zeros(1,numOfRows); % Binary flag for PDE detection.
allVarNames = [];
pdeVarNames = '';
dummys = 'DUMMYSPACE';
dummyt = 'DUMMYTIME';

if ( numOfRows == 1 ) % Not a system; can call convert2anon w/ two output args.

    [anFun indVarName allVarNames pdeVarNames eigVarNames commaSeparated] = ...
        setupLine(guifile, input{1}, type);

    % Check whether we have too many variables (this is a singl
    if ( ~isempty(pdeVarNames) ) % Not a PDE (or we can't do this type yet!)
        pdeflag = true;

        % Check whether we have a match of variable names, i.e. we want u
        % and u_t, not u and v_t:
        if ( (length(allVarNames) > 1) || (length(pdeVarNames) > 1) )
            error('Chebgui:setupFields:NumberOfVariables', ...
                'Too many variables.');
        end

        underScoreLocation = strfind(pdeVarNames{1},'_');
        if ( ~strcmp(allVarNames{1}, pdeVarNames{1}(1:underScoreLocation-1)) )
            error('Chebgui:setupFields:VariableNames', ...
                'Inconsistent variable names.');
        end
    end
    
    % Only need to convert the cell to a string for DE -- information has
    % been passed for BCs
    if ( strcmp(type, 'DE') )
        allVarString = allVarNames{1};
    end
    
    % Create the string which will become the anonymous function.
    % Put x (or t) as the first argument of the anonymous function if we
    % have a BVP or EIG.
    
    % !!! Check for a cell, make field a cell if we are in eigMode and
    % anFun is a cell.
    if ( strcmp(guifile.type, 'eig') && iscell(anFun) )
        field1 = ['@(', dummys, ',', allVarString, ') ', anFun{1}];
        field2 = ['@(', dummys, ',', allVarString, ') ', anFun{2}];
        field = {field1 ; field2};
    elseif ( ~strcmp(guifile.type, 'pde') && strcmp(type, 'DE') )
        field = ['@(', dummys, ',', allVarString, ') ', anFun];

    % Otherwise, add variables in front of what will be anonymous
    % functions. This is not needed for 'dirichlet','neumann',etc... This
    % can only happen in BCs, and in that case, allVarNames will be empty
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
    for k = 1:numOfRows
        [anFun{k} indVarNameTemp varNames pdeVarNames eigVarNames ...
            commaSeparated] = setupLine(guifile, input{k}, type);
        if ( strcmp(anFun{k}([1 end]), '[]') )
            anFun{k}([1 end]) = [];
        end

        allVarNames = [allVarNames ; varNames];

        % Only allow one time derivative in each line
        if ( length(pdeVarNames) > 1 )
            error('Chebgui:setupField:TooManyTimeDerivatives', ...
                'Only one time derivative per line allowed')
        end

        if ( isempty(pdeVarNames) )
            pdeVarNames = '|';
        else
            pdeflag(k) = 1;
        end
        
        if ( isempty(indVarName) )
            indVarName = indVarNameTemp;
        elseif ( ~mystrcmp(indVarName{1}, indVarNameTemp{1}) || ...
                ~mystrcmp(indVarName{2}, indVarNameTemp{2}) )
            error('Chebgui:Lexer:setupFields', ...
                'Different names for independent variables detected')
        end
        
        allPdeVarNames = [allPdeVarNames ; pdeVarNames];
    end

    % Remove duplicate variable names
    allVarNames = unique(allVarNames); 

    % Construct the handle part. For the DE field, we need to collect all
    % the variable names in one string. If we are working with BCs, we have
    % already passed that string in (as the parameter allVarString).
    if ( strcmp(type,'DE') )
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

        [ignored indx] = sort(indx); % Invert the index
        pdeflag = pdeflag(indx);     % Update the pdeflags
        allPdeVarNames(strcmp(allPdeVarNames, '|')) = []; % Delete the junk
        pdeVarNames = allPdeVarNames;
    end
    
    % If we are solving a BVP or EIG, we now need x as the first argument
    % as well. However, we don't want that variable in allVarString as we
    % use that information when setting up BCs. Create the string that goes
    % at the start of the final string.
    if ( (~strcmp(guifile.type, 'pde') && strcmp(type, 'DE')) || ...
            strcmp(type, 'BCnew') )
        fieldStart = ['@(', dummys, ',', allVarString, ') '];
    else
        fieldStart = ['@(', allVarString, ') '];
    end
        
    % Construct the function. Need to treat eig. problems separately as
    % there we can have nontrivial LHS and RHS at the same time.
    allAnFun = [];
    if ( strcmp(guifile.type, 'eig') && iscell(anFun{1}) )
        allAnFun1 = [];
        allAnFun2 = [];
        for k = 1:numOfRows
            allAnFun1 = [allAnFun1, anFun{indx(k)}{1}, ','];
            allAnFun2 = [allAnFun2, anFun{indx(k)}{2}, ','];
        end

        % Remove the last comma
        allAnFun1(end) = [];
        allAnFun2(end) = [];
        
        % Set up LHS and RHS fields
        field1 = [fieldStart, '[', allAnFun1, ']'];
        field2 = [fieldStart, '[', allAnFun2, ']'];
        field = {field1 ; field2};
    else
        for k = 1:numOfRows
            allAnFun = [allAnFun, anFun{indx(k)},  '; '];
        end
        allAnFun(end-1:end) = []; % Remove the last semicolon and space.
        
        field = [fieldStart, '[', allAnFun, ']'];
    end
    
end

end

function [field indVarName varNames pdeVarNames eigVarNames commaSeparated] ...
    = setupLine(guifile, input, type)

convertBCtoAnon = 0;

commaSeparated = 0; % Default value
% Create the variables x and t (corresponding to the linear function on the
% domain).
% x = xt; t = xt;
if ( ~isempty(strfind(input, '@')) ) % User supplied anon. function
    % Find the name of the independent variable (which we assume to be
    % either r, x or t)
    
    firstRPloc = strfind(input, ')');
    trimmedInput = input(firstRPloc+1:end);
    
    [field indVarName varNames pdeVarNames eigVarNames commaSeparated] = ...
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
   % Convert to anon. function string
    [field indVarName varNames pdeVarNames eigVarNames commaSeparated] = ...
        stringParser.str2anon(input, guifile.type, type);
end

end
    
% TODO:  Rename this function to something more sensible.
% Modified strcmp, if we compare with an empty string, give a match
function res =  mystrcmp(str1, str2)

if ( isempty(str1) || isempty(str2) )
    res = 1;
else
    res = strcmp(str1, str2);
end

end
