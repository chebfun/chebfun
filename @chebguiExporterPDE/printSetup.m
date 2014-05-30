function printSetup(fid, expInfo, guifile)

% Extract info from the expInfo struct
dom = expInfo.dom;
deInput = expInfo.deInput;

% PDE specific information
xName = expInfo.xName;
tName = expInfo.tName;
tt = expInfo.tt;
pdeflag = expInfo.pdeflag;
initInput = expInfo.initInput;

deString = expInfo.deString;
lbcInput = expInfo.lbcInput;
rbcInput = expInfo.rbcInput;

sol = expInfo.sol;
sol0 = expInfo.sol0;

allVarString = expInfo.allVarString;
allVarNames = expInfo.allVarNames;
periodic = expInfo.periodic;

fprintf(fid, '%% Create an interval of the space domain...\n');
fprintf(fid, 'dom = %s;\n',dom);
fprintf(fid, '%%...and a discretisation of the time domain:\n');
fprintf(fid, '%s = %s;\n',tName,tt);

fprintf(fid, '\n%% Make the right-hand side of the PDE.\n');
fprintf(fid, 'pdefun = %s;\n',deString);
if ( ~all(pdeflag) )
    fprintf(fid, ['pdeflag = [', num2str(pdeflag), ...
        ']; %% Zero when a variable is indep of time.\n']);
end

% Make assignments for left and right BCs.
fprintf(fid, '\n%% Assign boundary conditions.\n');
if ( ~isempty(lbcInput{1}) )
    lbcString = setupFields(guifile, lbcInput, 'BC', allVarString);
    fprintf(fid, 'bc.left = %s;\n', lbcString);
end

if ( ~isempty(rbcInput{1}) )
    rbcString = setupFields(guifile, rbcInput, 'BC', allVarString);
    fprintf(fid, 'bc.right = %s;\n', rbcString);
end

if ( periodic )
    fprintf(fid, 'bc = ''periodic'';\n');
end

% Set up the initial condition
fprintf(fid, '\n%% Construct a chebfun of the space variable on the domain,\n');
fprintf(fid, '%s = chebfun(@(%s) %s, dom);\n', xName, xName, xName);
if ( iscell(initInput) && (numel(initInput) > 1) )
    fprintf(fid, '%% and of the initial conditions.\n');
else
    fprintf(fid, '%% and of the initial condition.\n');
end

if ( (numel(deInput) == 1) && ~ischar(deInput) )
    % Get the strings of the dependant variable. Just use allVarNames.
    findx = strfind(initInput{1}, xName);
    initInput = vectorize(char(initInput));
    equalSign = find(initInput == '=', 1, 'last');
    if ( ~isempty(equalSign) )
        initInput = strtrim(initInput(equalSign+1:end));
    end
    if ( isempty(findx) )
        fprintf(fid, '%s = chebfun(%s,dom);\n', sol0, initInput);
    else
        fprintf(fid, '%s = %s;\n', sol0, vectorize(initInput));
    end
else
    % Get the strings of the dependant variables. Just use allVarNames
    s = allVarNames;
    
    % To deal with 'u = ...' etc in intial guesses
    order = [];
    guesses = [];
    inits = [];
    
    % Match LHS of = with variables in allVarNa
    for initCounter = 1:length(initInput)
        currStr = initInput{initCounter};
        equalSign = find(currStr == '=', 1, 'first');
        currVar = strtrim(currStr(1:equalSign-1));
        match = find(ismember(allVarNames, currVar) == 1);
        order = [order ; match];
        currInit = strtrim(currStr(1:equalSign-1));
        currGuess = vectorize(strtrim(currStr(equalSign+1:end)));
        guesses = [guesses ; {currGuess}];
        inits = [inits ; {currInit}];
    end
    [ignored, order] = sort(order);
    
    % If the initial guesses are all constants, we need to wrap them in a
    % chebfun call.
    for k = 1:numel(initInput)
        findx = strfind(initInput{k}, 'x');
        if ( ~isempty(findx) )
            break
        end
    end
    if ( isempty(findx) )
        for k = 1:numel(initInput)
            guesses{k} = sprintf('chebfun(%s,dom)', guesses{k});
        end
    end
    
    % These can be changed
    initText = '0';
    
    for k = 1:numel(initInput)
        fprintf(fid, '%s%s = %s;\n', inits{order(k)}, initText, ...
            guesses{order(k)});
    end
    fprintf(fid, '%s = [%s%s,', sol0, inits{order(1)}, initText);
    for k = 2:numel(initInput)-1
        fprintf(fid, ' %s%s,', inits{order(k)}, initText);
    end
    fprintf(fid, ' %s%s];\n', inits{order(end)}, initText);
    
end

end