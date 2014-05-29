function chebgui2mfile(exporter, guifile, pathname, filename)
% Export to .m file
%EXPORTBVP2MFILE    Export a BVP from CHEBGUI to a .m file.
%
%   See also: chebgui/export.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/chebfun/ for Chebfun information.

% TODO:  Documentation.


fullFileName = [pathname, filename];
fid = fopen(fullFileName, 'wt');

if ( ispc )
    userName = getenv('UserName');
else
    userName = getenv('USER');
end

fprintf(fid, ['%% %s - an executable M-file for solving a boundary value ' ...
    'problem.\n'], filename);
fprintf(fid, '%% Automatically created from chebfun/chebgui by user %s\n', ...
    userName);
fprintf(fid, '%% at %s on %s.\n\n', datestr(rem(now, 1), 13), ...
    datestr(floor(now)));

% Extract information from the GUI fields
dom = guifile.domain;
deInput = guifile.DE;
bcInput = guifile.BC;
initInput = guifile.init;

% Wrap all input strings in a cell (if they're not a cell already)
if ( isa(deInput, 'char') )
     deInput = cellstr(deInput);
end

if ( isa(bcInput, 'char') )
     bcInput = cellstr(bcInput);
end

if ( isa(initInput, 'char') )
     initInput = cellstr(initInput);
end

[deString, allVarString, indVarNameDE, dummy, dummy, dummy, allVarNames] = ...
    setupFields(guifile, deInput, 'DE');

% Do some error checking before we do further printing. Check that
% independent variable name match.
% Obtain the independent variable name appearing in the initial condition
useLatest = strcmpi(initInput{1}, 'Using latest solution');
if ( ~isempty(initInput{1}) && ~useLatest )
    [initString, dummy, indVarNameInit] = ...
        setupFields(guifile, initInput, 'BC', allVarString);
else
    indVarNameInit = {''};
end

% Make sure we don't have a discrepency in indVarNames
if ( ~isempty(indVarNameInit{1}) && ~isempty(indVarNameDE{1}) )
    if ( strcmp(indVarNameDE{1}, indVarNameInit{1}) )
        indVarNameSpace = indVarNameDE{1};
    else
        error('Chebgui:SolveGUIbvp', 'Independent variable names do not agree')
    end
elseif ( ~isempty(indVarNameInit{1}) && isempty(indVarNameDE{1}) )
    indVarNameSpace = indVarNameInit{1};
elseif ( isempty(indVarNameInit{1}) && ~isempty(indVarNameDE{1}) )
    indVarNameSpace = indVarNameDE{1};
else
    indVarNameSpace = 'x'; % Default value
end

% Replace the 'DUMMYSPACE' variable in the DE field
deString = strrep(deString, 'DUMMYSPACE', indVarNameSpace);
deString = prettyprintfevalstring(deString, allVarNames);

% Support for periodic boundary conditions
if ( ~isempty(bcInput{1}) && strcmpi(bcInput{1}, 'periodic') )
    bcInput{1} = [];
    periodic = true;
else
    periodic = false;
end

% Print the BVP
fprintf(fid, '%% Solving\n');
for k = 1:numel(deInput)
    fprintf(fid, '%%   %s,\n', deInput{k});
end
fprintf(fid, '%% for %s in %s', indVarNameSpace, dom);
if ( ~isempty(bcInput{1}) )
    fprintf(fid, ', subject to\n%%');
    for k = 1:numel(bcInput)
        fprintf(fid, '   %s', bcInput{k});
        if ( (k ~= numel(bcInput)) && (numel(bcInput) > 1) )
            fprintf(fid, ',\n%%');
        end
    end
    fprintf(fid, '.\n');
elseif periodic
    fprintf(fid, ', subject to periodic boundary conditions.\n\n');
else
    fprintf(fid, '.\n');
end

fprintf(fid, '\n%% Define the domain.\n');
fprintf(fid, 'dom = %s;\n', dom);
fprintf(fid, ['\n%% Assign the differential equation to a chebop on that ' ...
    'domain.\n']);
fprintf(fid, 'N = chebop(%s,dom);\n', deString);

% Setup for the rhs
fprintf(fid, ['\n%% Set up the rhs of the differential equation so that ' ...
    'N(%s) = rhs.\n'], allVarString);

% If we have a coupled system, we need create a array of the inputs
if ( size(deInput, 1) > 1 )
    deRHSprint = ['['];
    for counter = 1:size(deInput,1)
        deRHSprint = [deRHSprint num2str(0) ','];
    end
    deRHSprint(end) = []; % Remove the last comma
    deRHSprint = [deRHSprint, ']'];
else
    deRHSprint = num2str(0);
end
fprintf(fid, 'rhs = %s;\n', deRHSprint);

% Make assignments for BCs.
fprintf(fid, '\n%% Assign boundary conditions to the chebop.\n');
if ( ~isempty(bcInput{1}) )
    bcString = setupFields(guifile, bcInput, 'BCnew', allVarString );
    bcString = strrep(bcString, 'DUMMYSPACE', indVarNameSpace);
    bcString = prettyprintfevalstring(bcString, allVarNames);
    fprintf(fid, 'N.bc = %s;\n', bcString);
end
if ( periodic )
    fprintf(fid, 'N.bc = ''periodic'';\n');
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
%     fprintf(fid,'N.init = %s;\n',vectorize(char(initInput)));
    initInput = cellstr(initInput);
    if ( numel(initInput) == 1 )
        guessInput = vectorize(strtrim(char(initInput{1})));
        equalSign = find(guessInput == '=', 1, 'last');
        if ( ~isempty(equalSign) )
            guessInput = guessInput(equalSign+1:end); 
        end
        fprintf(fid, 'N.init = %s;\n', guessInput);
    else
        % To deal with 'u = ...' etc in intial guesses
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

        [ignored, order] = sort(order);
        initText = '_init';
        for k = 1:numel(initInput)
            fprintf(fid, '%s%s = %s;\n', inits{order(k)}, initText, ...
                guesses{order(k)});
        end
        fprintf(fid, 'N.init = [%s%s,', inits{order(1)}, initText);
        for k = 2:numel(initInput)-1
            fprintf(fid, ' %s%s,', inits{order(k)}, initText);
        end
        fprintf(fid, ' %s%s];\n', inits{order(end)}, initText);

    end
end

% Set up preferences
fprintf(fid, '\n%% Setup preferences for solving the problem.\n');
fprintf(fid, 'options = cheboppref();\n');

% Option for tolerance
tolInput = guifile.tol;

if ( ~isempty(tolInput) )
    fprintf(fid, '\n%% Option for tolerance.\n');
    fprintf(fid, 'options.errTol = %s;\n', tolInput);
%     fprintf(fid,'options.restol = %s;\n',tolInput);
end

fprintf(fid, 'options.display = ''iter'';\n');

% Option for damping
dampedOnInput = guifile.options.damping;

fprintf(fid, '\n%% Option for damping.\n');
if ( strcmp(dampedOnInput, '1') )
    fprintf(fid, 'options.damped = ''on'';\n');
else
    fprintf(fid, 'options.damped = ''off'';\n');
end

fprintf(fid, '\n%% Option for discretization (either @colloc2 or @ultraS).\n');
fprintf(fid, 'options.discretization = @%s;\n', ...
    func2str(guifile.options.discretization));

% Option for plotting
plottingOnInput = guifile.options.plotting;

if ( ~strcmp(plottingOnInput, 'off') )
    fprintf(fid, ['\n%% Option for determining how long each Newton step ' ...
        'is shown.\n']);
    fprintf(fid, 'options.plotting = %s;\n', plottingOnInput);
else
    fprintf(fid, ['\n%% Option for determining how long each Newton step ' ...
        'is shown.\n']);
    fprintf(fid,'options.plotting = ''off'';\n');
end

fprintf(fid, ['\n%% Solve the problem using solvebvp (a routine which ' ...
    'offers the same\n%% functionality as nonlinear backslash, but with '...
    'more customizability).\n']);
% fprintf(fid,'[u normVec] = solvebvp(N,rhs,''options'',options);\n');
fprintf(fid, 'u = solvebvp(N,rhs,options);\n');

fprintf(fid, '\n%% Create a plot of the solution.\n');

fprintf(fid, ['figure\nplot(u,''LineWidth'',2)\ntitle(''Final solution'')', ...
    ', xlabel(''%s'')'], indVarNameSpace);
if ( numel(allVarNames) == 1 )
    fprintf(fid, ', ylabel(''%s'')', allVarNames{:});
else
    leg = '';
    for k = 1:numel(allVarNames)-1
        leg = [leg '''' allVarNames{k} '''' ','];
    end
    leg = [leg '''' allVarNames{k+1} ''''];
    fprintf(fid, ', legend(%s)\n', leg);
end

fclose(fid);

end


function str = prettyprintfevalstring(str, varnames)

for k = 1:numel(varnames)
    oldstr = ['feval(' varnames{k} ','];
    newstr = [varnames{k} '('];
    str = strrep(str, oldstr, newstr);
    oldstr = [varnames{k} '(''end'''];
    newstr = [varnames{k} '(end'];
    str = strrep(str, oldstr, newstr);
    oldstr = [varnames{k} '(''right'''];
    newstr = [varnames{k} '(end'];
    str = strrep(str, oldstr, newstr);
    oldstr = [varnames{k} '(''start'''];
    newstr = [varnames{k} '(' varnames{k} '.ends(1)'];
    str = strrep(str, oldstr, newstr);
    oldstr = [varnames{k} '(''left'''];
    newstr = [varnames{k} '(' varnames{k} '.ends(1)'];
    str = strrep(str, oldstr, newstr);
end

end