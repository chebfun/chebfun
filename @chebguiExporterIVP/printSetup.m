function printSetup(fid, expInfo, guifile)
%PRINTSETUP   Print commands for setting up problem.
%
% Calling sequence:
%   PRINTPOSTSOLVER(FID, EXPINFO)
% where
%   FID:        ID of a file-writing stream.
%   EXPINFO:    Struct containing information for printing the problem.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract info from the EXPINFO struct:
dom = expInfo.dom;
deString = expInfo.deString;
deInput = expInfo.deInput;
bcInput = expInfo.bcInput;
initInput = expInfo.initInput;
allVarString = expInfo.allVarString;
allVarNames = expInfo.allVarNames;
indVarNameSpace = expInfo.indVarNameSpace;
useLatest = expInfo.useLatest;
timeSteppingSolver = expInfo.timeSteppingSolver;

% Print commands for problem set-up:
fprintf(fid, '\n%%%% Problem set-up.');
fprintf(fid, '\n%% Define the domain.\n');
fprintf(fid, 'dom = %s;\n', dom);
fprintf(fid, ['\n%% Assign the differential equation to a chebop on that ' ...
    'domain.\n']);
fprintf(fid, 'N = chebop(%s, dom);\n', deString);

% Setup for the RHS
fprintf(fid, ['\n%% Set up the rhs of the differential equation so that ' ...
    'N(%s) = rhs.\n'], allVarString);

% If we have a coupled system, we need create a array of the inputs:
if ( size(deInput, 1) > 1 )
    deRHSprint = '[';
    for counter = 1:size(deInput,1)
        deRHSprint = [deRHSprint num2str(0) ',']; %#ok<AGROW>
    end
    deRHSprint(end) = []; % Remove the last comma
    deRHSprint = [deRHSprint, ']'];
    
else
    deRHSprint = num2str(0);
    
end
fprintf(fid, 'rhs = %s;\n', deRHSprint);

% Make assignments for BCs.

% Are we dealing with an initial or final value problem?
isIorF = isIVPorFVP(guifile, allVarNames);

% Since we're solving an IVP or an FVP, we need to reformulate the BC input, so
% that it can be assigned to the N.LBC or N.RBC fields:
newBCinput = chebgui.bcReform(expInfo.dom, bcInput, isIorF);
bcString = setupFields(guifile, newBCinput, 'BC', allVarString);

% Make assignments for BCs.
fprintf(fid, '\n%% Assign boundary conditions to the chebop.\n');
if ( isIorF == 1 )
    fprintf(fid, 'N.lbc = %s;\n', bcString);
else
    fprintf(fid, 'N.rbc = %s;\n', bcString);    
end

% Set up for the initial guess of the solution.
if ( useLatest )
    fprintf(fid, ['\n%% Note that it is not possible to use the "Use ' ...
        'latest" option \n%% when exporting to .m files. \n']);
    
elseif ( ~isempty(initInput{1}) )
    fprintf(fid, '\n%% Construct a linear chebfun on the domain, \n');
    fprintf(fid, '%s = chebfun(@(%s) %s, dom);\n', ...
        indVarNameSpace, indVarNameSpace, indVarNameSpace);
    fprintf(fid, '%% and assign an initial guess to the chebop.\n');
    initInput = cellstr(initInput);
    
    if ( numel(initInput) == 1 )
        % Only one input passed (scalar problem).
        guessInput = vectorize(strtrim(char(initInput{1})));
        equalSign = find(guessInput == '=', 1, 'last');
        if ( ~isempty(equalSign) )
            guessInput = guessInput(equalSign+1:end);
        end
        fprintf(fid, 'N.init = %s;\n', guessInput);
        
    else
        % To deal with 'u = ...' etc in initial guesses
        order = [];
        guesses = [];
        inits = [];
        
        % Match LHS of = with variables in allVarName
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
        
        % Sort the variables to ensure they get assigned the correct guesses.
        [ignored, order] = sort(order);
        initText = '_init';
        for k = 1:numel(initInput)
            fprintf(fid, '%s%s = %s;\n', inits{order(k)}, initText, ...
                guesses{order(k)});
        end
        fprintf(fid, 'N.init = [%s%s;', inits{order(1)}, initText);
        for k = 2:numel(initInput)-1
            fprintf(fid, ' %s%s;', inits{order(k)}, initText);
        end
        fprintf(fid, ' %s%s];\n', inits{order(end)}, initText);
    end
    
end
