function varargout = solveGUIivp(guifile, handles)
%SOLVEGUIIVP     Solve an IVP, specified by a CHEBGUI object.
%
% Calling sequence:
%
%   VARARGOUT = SOLVEGUIIVP(GUIFILE, HANDLES)
%
% where
%   
%   GUIFILE:    A CHEBGUI object, describing the problem.
%   HANDLES:    A MATLAB handle to the chebguiwindow figure.
%
% If the method is called by pressing the 'Solve' button on the GUI,
%   VARARGOUT{1}:   Will be a MATLAB handle to the chebguiwindow figure, which
%                   has been updated to contain the solution and other useful
%                   results for the problem solved.
%
% If the method is called by calling the command explicitly with a CHEBGUI
% object (e.g. [U, INFO] = SOLVEGUIVP(GUIFILE) from the command line),
%   VARARGOUT{1}:   The solution to the problem specified by GUIFILE.
%   VARARGOUT{2}:   The INFO struct returned by the chebop/solveivp() method.
%
% See also: chebgui/solveGUI, chebgui/solveGUIbvp.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Handles will be an empty variable if we are solving without using the GUI
if ( nargin < 2 )
    guiMode = 0;
else
    guiMode = 1;
end

% Call the exportInfo method of the chebguiExporterIVP class, which takes care
% of extracting most information we need from GUIFILE.
expInfo = chebguiExporterIVP.exportInfo(guifile);

dom = str2num(expInfo.dom);
deString = expInfo.deString;
allVarString = expInfo.allVarString;
allVarNames = expInfo.allVarNames;
indVarNameSpace = expInfo.indVarNameSpace;
bcInput = expInfo.bcInput;
initInput = expInfo.initInput;
% Create the independent variable on DOM.
xt = chebfun(@(x) x, dom);

% Assign the variables names to HANDLES.
handles.varnames = allVarNames;
handles.indVarName = {indVarNameSpace};

% If we only have one variable appearing in allVarNames, the problem is a
% scalar problem.
scalarProblem = length(allVarNames) == 1;

% Are we dealing with an initial, final, or boundary value problem?
isIorF = isIVPorFVP(guifile, expInfo.allVarNames);

% Assert that we're actually solving an IVP!
assert(logical(isIorF), 'CHEBFUN:CHEBGUI:solveGUIivp:problemType', ...
    ['Problem detected to be a boundary-value problem, but attempted to ' ...
    'solve as an initial value problem. Please check problem specification']);
    
% Reformulate the boundary condition so that we can assign it to the .LBC or
% .RBC field below.
newBCinput = chebgui.bcReform(expInfo.dom, bcInput, isIorF);
bcString = setupFields(guifile, newBCinput, 'BC', allVarString);

% Set up the initial guesses (which might be useful if for some reason, we are
% asking to solve the IVP globally with Newton's method).
guess = [];

% Find whether the user wants to use the latest solution as a guess. This is
% only possible when calling from the GUI
if ( guiMode )
    useLatest = strcmpi(get(handles.input_GUESS, 'String'), ...
        'Using latest solution');
    if ( useLatest )
        guess = handles.latest.solution;
    end
end

% Assign r, x or t as the linear function on the domain if indVarName is
% not empty
eval([indVarNameSpace, '=xt;']);

% Convert the string to proper anon. function using eval
DE = eval(deString);

% Do the same in the .BC field if it isn't empty
if ( ~isempty(bcString) )
    bcString = strrep(bcString, 'DUMMYSPACE', indVarNameSpace);
    BC = eval(bcString);
end

% Convert the initial condition passed to a CHEBFUN.
if ( ~isempty(initInput{1}) && isempty(guess) )
    if ( iscellstr(initInput) )
        order = [];
        guesses = [];

        % Match LHS of = with variables in allVarNames
        for initCounter = 1:length(initInput)
            currStr = initInput{initCounter};
            equalSign = find(currStr == '=');
            currVar = strtrim(currStr(1:equalSign-1));
            match = find(ismember(allVarNames, currVar) == 1);
            order = [order ; match];
            currGuess = strtrim(currStr(equalSign+1:end));
            guesses = [guesses ; {currGuess}];
        end
        
        % If order is still empty, that means that initial guess were
        % passed on the form '2*x', rather than 'u = 2*x'. Allow that for
        % scalar problems, throw an error otherwise.
        if ( isempty(order) && scalarProblem )
            guess = eval(vectorize(initInput{1}));
            % If we have a scalar guess, convert to a chebfun
            if ( isnumeric(guess) )
                guess = 0*xt + guess;
            end
        elseif ( length(order) == length(guesses) )
            % We have a guess to match every dependent variable in the
            % problem.
            guess = cell(length(order),1);
            for guessCounter = 1:length(guesses)
                guessLoc = find(order == guessCounter);
                tempGuess = eval(vectorize(guesses{guessLoc}));
                if ( isnumeric(tempGuess) )
                    tempGuess = 0*xt + tempGuess;
                end
                guess{guessCounter} = tempGuess;
            end
            
            guess = chebmatrix(guess);
        else % throw an error
            error('CHEBFUN:CHEBGUI:solveGUIivp:initialGuess', ...
                ['Error constructing initial guess.  Please make sure ' ...
                 'guesses are of the form u = 2*x, v = sin(x), ...']);
        end
    else
        % initInput is not a cell string, must have only received one line.
        guessInput = vectorize(initInput);
        equalSign = find(guessInput == '=');
        if ( isempty(equalSign) )
            equalSign = 0;
        end
        guessInput = guessInput(equalSign+1:end);
        guess =  chebfun(guessInput, [a b]);
    end
end

% Before continuing, clear any previous information about the norm of the Newton
% updates in the handles struct:
handles.normDelta = [];

% Create the CHEBOP.
N = chebop(DE, dom);

% Assign boundary conditions, depending on whether we're solving an initial
% value problem or a final value problem:
if ( isIorF == 1 )
    N.lbc = BC;
else
    N.rbc = BC;
end

% Assign initial guess if it was passed
if ( ~isempty(guess) )
    N.init = guess;
end

% Obtain the CHEBOPPREF to pass to the solvers.
options = setupODEoptions(handles.guifile, expInfo);

% Are we solving the problem globally, or with one of the MATLAB solvers?
solvingGlobally = isempty(strfind(func2str(options.ivpSolver),'chebfun'));
% Various things we only need to think about when in the GUI, changes GUI
% compenents.
if ( guiMode )
    set(handles.iter_list, 'String', '');
    set(handles.iter_text, 'Visible', 'On');
    set(handles.iter_list, 'Visible', 'On');
    
    xLimit = dom([1 end]);
    handles.xLim = xLimit;

    set(handles.fig_sol, 'Visible', 'On');
    set(handles.fig_norm, 'Visible', 'On');
    
    set(handles.popupmenu_bottomFig, 'Value', 1);
end

% Call solveivp with different arguments depending on whether we're in GUI
% or not. If we're not in GUI mode, we can finish here.
if ( guiMode )
    if ( solvingGlobally )        
        % We use the BVP displayInfo method, as the only time it will be used is
        % if we're solving the problem globally, in which case, it's the
        % appropriate display method.
        displayFunction = ...
            @(mode, varargin) chebgui.displayBVPinfo(handles, mode, varargin{:});
        
        [u, info] = solvebvp(N, 0, options, displayFunction);
    else
        % Solving using one of the MATLAB solvers.
        u = solveivp(N, 0, options);
        info = struct('isLinear', 1, 'error', 0);
        if (isa(u,'chebfun'))
            len = length(u);
        else
            len = max(cellfun(@length, u.blocks));
        end
        
        if ( isIorF == 1 )
            detStr = 'Initial value problem detected';
        else
            detStr = 'Final value problem detected';
        end
        
        str = {detStr, ; 
            sprintf('Number of pieces of the solution: %i.', ...
            length(u.domain)-1);
            sprintf('Total length of solution: %i.', len)};
        set(handles.iter_list, 'String',  str)
        set(handles.iter_list, 'Value', 1)
    end
else
    % Not in GUI mode.
    
    [u, info] = mldivide(N, 0, options);
    varargout{1} = u;
    varargout{2} = info;
end

% ISLINEAR is returned as a vector in the INFO structure (with elements
% corresponding to DE and BCs. Convert to a binary, 1 if everything in the
% problem is linear, 0 otherwise
isLinear = all(info.isLinear);

% Now do some post solving processing, specific to running in GUI mode.
if ( guiMode )
    % Convert u from a CHEBMATRIX to an array valued CHEBFUN to make plotting
    % easier.
    if ( isa(u, 'chebmatrix') )
        u = chebfun(u);
    end
    
    % Store in handles latest chebop, solution, vector of norm of updates etc.
    % (enables exporting later on)
    handles.latest.type = 'ivp';
    handles.latest.solution = u;
    handles.latest.norms = info.error;
    handles.latest.chebop = N;
    handles.latest.options = options;
    
    % Notify the GUI we have a solution available
    handles.hasSolution = 1;
    
    % Plot
    axes(handles.fig_sol)
    plot(u, 'Linewidth',2)
    % Do different things depending on whether the solution is real or not
    if ( isreal(u) )
        xlim(xLimit)
    else
        axis equal
    end

    % Only show legend if we were solving a coupled system:
    if ( length(allVarNames) > 1 )
        legend(allVarNames)
    end

    % Show grid?
    if ( guifile.options.grid )
        grid on
    end

    % Different titles of top plot if we had a linear problem:
    if ( solvingGlobally && ~isLinear )
        title('Solution at end of iteration');
    else
        title('Solution');
    end

    % If we were solving a nonlinear problem, we show a plot of the norm of the
    % Newton updates after solution has been found. For a linear problem, we
    % show a PLOTCOEFFS plot.
    if ( solvingGlobally && ~isLinear )
        % Store the norm of the Newton updates
        handles.normDelta = info.normDelta;

        axes(handles.fig_norm)
        normDelta = info.normDelta;
        semilogy(normDelta, '-*', 'Linewidth', 2)
        title('Norm of updates')
        xlabel('Iteration number')

        if ( length(normDelta) > 1 )
            XTickVec = 1:max(floor(length(normDelta)/5), 1):length(normDelta);
            set(gca, 'XTick', XTickVec)
            xlim([1 length(normDelta)])
            grid on
        else % Don't display fractions on iteration plots
            set(gca,'XTick', 1)
        end
    else
        % If we're solving a linear problem, or using the MATLAB solvers, plot
        % the Chebyshev coefficients instead.
        axes(handles.fig_norm)
        plotcoeffs(u, 'linewidth', 2)
        title('Chebyshev coefficients of the solution')
        set(handles.popupmenu_bottomFig, 'Value', 2);
        grid on
    end
    
    % Return the handles as varargout.
    varargout{1} = handles;
end

end
