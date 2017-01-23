function varargout = solveGUIbvp(guifile, handles)
%SOLVEGUIBVP   Solve a BVP, specified by a CHEBGUI object.
%
% Calling sequence:
%
%   VARARGOUT = SOLVEGUIBVP(GUIFILE, HANDLES)
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
% object (e.g. [U, INFO] = SOLVEGUIBVP(GUIFILE) from the command line),
%   VARARGOUT{1}:   The solution to the problem specified by GUIFILE.
%   VARARGOUT{2}:   The INFO struct returned by the chebop/solvebvp method.
%
% See also: chebgui/solveGUI, chebgui/solveGUIivp.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Handles will be an empty variable if we are solving without using the GUI
if ( nargin < 2 )
    guiMode = 0;
else
    guiMode = 1;
end

% Call the exportInfo method of the chebguiExporterBVP class, which takes care
% of extracting most information we need from GUIFILE.
expInfo = chebguiExporterBVP.exportInfo(guifile);

dom = str2num(expInfo.dom);
deString = expInfo.deString;
allVarString = expInfo.allVarString;
allVarNames = expInfo.allVarNames;
indVarNameSpace = expInfo.indVarNameSpace;
bcInput = expInfo.bcInput;
periodic = expInfo.periodic;
initInput = expInfo.initInput;
numVars = expInfo.numVars;
% Create the independent variable on DOM.
xt = chebfun(@(x) x, dom);

% Assign the variables names to HANDLES.
handles.varnames = allVarNames;
handles.indVarName = {indVarNameSpace};

% Obtain the boundary conditions to be imposed.
if ( periodic )
    bcString = '';
    BC = 'periodic';
elseif ( isempty(bcInput{1}) )
    bcString = '';
    BC = [];
else
    bcString = setupFields(guifile, bcInput, 'BCnew', allVarString);
end

% Set up the initial guesses.
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

% If we got an initial guess passed, convert it to a CHEBFUN:
if ( ~isempty(initInput{1}) && isempty(guess) )
    guess = chebgui.constructInit(initInput, allVarNames, indVarNameSpace, xt);
end

% Before continuing, clear any previous information about the norm of the Newton
% updates in the handles struct:
handles.normDelta = [];

% Create the CHEBOP, including boundary conditions:
N = chebop(DE, dom, BC);

% Assign initial guess if it was passed
if ( ~isempty(guess) )
    N.init = guess;
end

% Obtain the CHEBOPPREF to pass to the solvers.
options = setupODEoptions(handles.guifile, expInfo);


% Various things we only need to think about when in the GUI, changes GUI compenents.
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

% Call solvebvp with different arguments depending on whether we're in GUI
% or not. If we're not in GUI mode, we can finish here.
if ( guiMode )
    displayFunction = ...
        @(mode, varargin) chebgui.displayBVPinfo(handles, mode, varargin{:});
    [u{1:numVars,1}, info] = mldivide(N, 0, options, displayFunction);
else
    [u{1:numVars,1}, info] = mldivide(N, 0, options);
    varargout{1} = u;
    varargout{2} = info;
end

% Convert the cell-array we get returned above to a chebmatrix
u = chebmatrix(u);

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
    handles.latest.type = 'bvp';
    handles.latest.solution = u;
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
    if ( ~isLinear )
        set(handles.panel_figSol, 'title', 'Solution at end of iteration')
    else
        set(handles.panel_figSol, 'title', 'Solution')
    end

    % If we were solving a nonlinear problem, we show a plot of the norm of the
    % Newton updates after solution has been found. For a linear problem, we
    % show a PLOTCOEFFS plot.
    if ( ~isLinear )
        % Store the norm of the Newton updates
        handles.latest.normDelta = info.normDelta;

        axes(handles.fig_norm)
        normDelta = info.normDelta;
        semilogy(normDelta, '-*', 'Linewidth', 2)
        set(handles.panel_figNorm, 'title', 'Norm of updates')
        set(handles.fig_norm, 'fontsize', handles.fontsizePanels);

        if ( length(normDelta) > 1 )
            XTickVec = 1:max(floor(length(normDelta)/5), 1):length(normDelta);
            set(gca, 'XTick', XTickVec)
            xlim([1 length(normDelta)])
            grid on
        else % Don't display fractions on iteration plots
            set(gca,'XTick', 1)
        end
    else
        axes(handles.fig_norm)
        % Show the plotcoeffs plot. Grab its title, set as the title of the
        % panel, then hide the title of the plot and the x/ylabels (to avoid
        % issues at large fontsizes):
        plotcoeffs(u, 'linewidth', 2)
        plotCoeffsTitle = get(get(handles.fig_norm, 'title'), 'String');
        set(handles.panel_figNorm, 'title', plotCoeffsTitle);
        title('');
        xlabel('');
        ylabel('');
        
        set(handles.fig_norm, 'fontsize', handles.fontsizePanels);
        % Store an empty vector for the norm of the updates (since problem was
        % linear)
        handles.latest.normDelta = [];
        % In older versions of MATLAB, need to change the title font-size manually:
        if verLessThan('matlab', '8.4')
            set(get(handles.fig_norm, 'title'), 'fontsize', handles.fontsizePanels)
        else
            set(get(handles.fig_norm, 'title'), 'fontweight', 'normal')
        end
        set(handles.popupmenu_bottomFig, 'Value', 2);
        grid on
    end

    % Update the fontsize of plots
    set(handles.fig_sol, 'fontsize', handles.fontsizePanels);
    set(handles.fig_norm, 'fontsize', handles.fontsizePanels);

    % Return the handles as varargout.
    varargout{1} = handles;
end

end
