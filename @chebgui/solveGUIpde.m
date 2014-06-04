function varargout = solveGUIpde(guifile, handles)
%SOLVEGUIPDE   Solve a PDE, specified by a CHEBGUI object.
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
%   VARARGOUT{1}:   The time range of the problem specified by GUIFILE.
%   VARARGOUT{2}:   The a CHEBMATRIX containing the solution returned by pde15s.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Handles will be an empty variable if we are solving without using the GUI
if ( nargin < 2 )
    guiMode = 0;
else
    guiMode = 1;
end

% Extract information from the guifile
deInput = guifile.DE;
lbcInput = guifile.LBC;
rbcInput = guifile.RBC;
initInput = guifile.init;
dom = str2num(guifile.domain);
if ( length(dom) > 2 )
    warning('Chebgui:SolveGUIpde', ...
        ['PDE solver does not accept domains with breakpoints. ' ...
         'Breakpoints are being ignored.']);
    dom = dom([1 end]);
end
tolInput = guifile.tol;
tt = str2num(guifile.timedomain);

% Wrap all input strings in a cell (if they're not a cell already)
if ( isa(deInput,'char') )
    deInput = cellstr(deInput);
end

if ( isa(lbcInput,'char') )
    lbcInput = cellstr(lbcInput);
end

if ( isa(rbcInput,'char') )
    rbcInput = cellstr(rbcInput);
end

if ( isa(initInput,'char') )
    initInput = cellstr(initInput);
end

% Convert the input to the an. func. format, get information about the
% independent variable in the problem.
[deString, allVarString, indVarNameDE, dummy, pdeflag, dummy, allVarNames] ...
    = setupFields(guifile, deInput, 'DE');
handles.varnames = allVarNames;
if ( ~any(pdeflag) )
    s = ['Input does not appear to be a PDE, or at least is not a ' ...
         'supported type. Perhaps you need to switch to ''ODE'' mode?'];
    error('Chebgui:SolveGUIpde', s);
end

% Are we solving a scalar problem?
scalarProblem = length(allVarNames) == 1;

% Obtain the independent variable name appearing in the initial condition
if ( ~isempty(initInput{1}) )
    % If we have a scalar problem, we're OK with no dependent variables
    % appearing in the initial guess
    if ( scalarProblem ) 
        [dummy, dummy, indVarNameInit] = ...
            setupFields(guifile, initInput, 'INITSCALAR', allVarString);
        % If the initial guess was just passed a constant, indVarNameInit will
        % be empty. For consistency (allowing us to do the if-statement below),
        % convert to an empty cell.
        if ( isempty(indVarNameInit) )
            indVarNameInit = {''};
        end       
    else
        [dummy, dummy, indVarNameInit] = ...
            setupFields(guifile, initInput, 'INIT', allVarString);        
    end
else
    indVarNameInit = {''};
end

% Make sure we don't have a disrepency in indVarNames. Create a new variable
% indVarName which contains both independent variable, the first entry
% corresponds to space, the second to time.
if ( ~isempty(indVarNameInit{1}) && ~isempty(indVarNameDE{1}) )
    if ( strcmp(indVarNameDE{1}, indVarNameInit{1}) )
        indVarName{1} = indVarNameDE{1};
    else
        error('Chebgui:SolveGUIpde', 'Independent variable names do not agree')
    end
elseif ( ~isempty(indVarNameInit{1}) && isempty(indVarNameDE{1}) )
    indVarName{1} = indVarNameInit{1};
elseif ( isempty(indVarNameInit{1}) && ~isempty(indVarNameDE{1}) )
    indVarName{1} = indVarNameDE{1};
else
    indVarName{1} = 'x'; % Default value
end

if ( ~isempty(indVarNameDE{2}) )
    % Find what the name of the time variable specified is:
    indVarName{2} = indVarNameDE{2};
else
    error('Chebgui:SolveGUIpde', ...
        'No time variable specified, please do so via using subscript.');
end

% Did we use the same variable name for both space and time?
if ( strcmp(indVarName{1}, indVarName{2}) )
     error('Chebgui:SolveGUIpde', ...
        'The same variable appears to be used as space and time variable');
end
handles.indVarName = indVarName;

% Create a string with the variables used in the problem
variableString = [indVarName{2}, ',', indVarName{1}, ','];

% String that will be converted to an anonymous function:
deString = ['@(' variableString, deString(3:end)];

% Convert the string to proper anon. function using eval
DE = eval(deString);

% Support for periodic boundary conditions
if ( (~isempty(lbcInput{1}) && strcmpi(lbcInput{1}, 'periodic')) || ...
        (~isempty(rbcInput{1}) && strcmpi(rbcInput{1}, 'periodic')) )
    lbcInput{1} = [];
    rbcInput{1} = [];
    periodic = true;
else
    periodic = false;
end

% Do we have a left BC specified?
if ( ~isempty(lbcInput{1}) )
    lbcString = setupFields(guifile, lbcInput, 'BC', allVarString);
    
    idx = strfind(lbcString, ')');
    if ( ~isempty(idx) )
        
        % Inject the t and x variables into the anonymous function:
        lbcString = [lbcString(1:2), variableString, lbcString(3:end)];

    end

    LBC = eval(lbcString);
else
    LBC = [];
end

% Do we have a right BC specified?
if ( ~isempty(rbcInput{1}) )
    rbcString = setupFields(guifile, rbcInput, 'BC', allVarString);

    idx = strfind(rbcString, ')');
    if ( ~isempty(idx) )
        
        % Inject the t and x variables into the anonymous function:
        rbcString = [rbcString(1:2), variableString, rbcString(3:end)];
        
    end
    
    RBC = eval(rbcString);
else
    RBC = [];
end

% No boundary conditions are a no-go:
if ( isempty(lbcInput) && isempty(rbcInput) )
    error('chebfun:bvpgui','No boundary conditions specified');
end

% Boundary condition
if ( periodic )
    bc = 'periodic';
else
    bc = [];
    bc.left = LBC;
    bc.right = RBC;
end

% Start gathering options. Start by obtaining a PDE options struct.
opts = pdeset;
defaultTol = opts.Eps;
defaultLineWidth = 2;

% Check that we're not imposing too strict tolerance:
if ( isempty(tolInput) )
    tolNum = defaultTol;
else
    tolNum = str2num(tolInput);
end
% Find out what the current CHEBFUN eps is set at.
cpref = chebfunpref();
if ( tolNum < cpref.eps )
    warndlg('Tolerance specified is less than current chebfun epsilon', ...
        'Warning','modal');
    uiwait(gcf)
end

tol = tolNum;

% Set up the initial condition
if ( iscellstr(initInput) )
    % We have multiple lines specifying the initial condition:
    order = [];
    inits = [];

    % Match LHS of = with variables in allVarName
    for initCounter = 1:length(initInput)
        currStr = initInput{initCounter};
        equalSign = find(currStr == '=');
        currVar = strtrim(currStr(1:equalSign-1));
        match = find(ismember(allVarNames, currVar) == 1);
        order = [order ; match];
        currInit = strtrim(currStr(equalSign+1:end));
        inits = [inits ; {currInit}];
    end
    
    u0 = chebfun;
    
    if ( isempty(order) && scalarProblem )
        u0 = chebfun(vectorize(initInput{1}), dom);
        u0 = simplify(u0, tol);
    else
        for k = 1:length(inits)
            initLoc = find(order == k);
            init_k = vectorize(inits{initLoc});
            u0k = chebfun(init_k, dom);
            u0k = simplify(u0k, tol);
            u0 = [u0, u0k];
        end
    end
else
    % Only one line specifiying the initial condition (for a scalar problem):
    initInput = vectorize(initInput);
    
    % Is the input of the form '3' or 'u=3'?
    equalSign = find(initInput == '=');
    if ( isempty(equalSign) )
        equalSign = 0;
    end
    initInput = initInput(equalSign+1:end);
    u0 =  chebfun(initInput, dom);
    u0 = simplify(u0, tol);
end


% Get the GUI ready for plotting.
if ( guiMode )
    set(handles.fig_sol, 'Visible', 'On');
    set(handles.fig_norm, 'Visible', 'On');
    cla(handles.fig_sol, 'reset')
    cla(handles.fig_norm, 'reset')
    handles.gui = 1;
else
    handles.gui = 0;
end

% Continue setting up options.

% Give the OPTS structure the handle to the chebguiwindow figure, so that pde15s
% can plot to it:
opts.handles = handles;

% Be default, don't do a hold on plot.
opts.HoldPlot = false;
opts.Eps = tolNum;
if ( ~all(pdeflag) )
    opts.PDEflag = pdeflag;
end

opts.Plot = guifile.options.plotting;

if ( guifile.options.pdeholdplot )
    opts.HoldPlot = 'on';
else
    opts.HoldPlot = 'off';
end

% Do we want to fix the y-limits while solving the PDE?
ylim1 = guifile.options.fixYaxisLower;
ylim2 = guifile.options.fixYaxisUpper;
if ( ~isempty(ylim1) && ~isempty(ylim2) )
    opts.YLim = [str2double(ylim1), str2double(ylim2)];
end
opts.PlotStyle = {'LineWidth', defaultLineWidth};

% Options for fixed N
if ( ~isempty(guifile.options.fixN) )
    opts.N = str2num(guifile.options.fixN);
end

% Try solving the PDE!
try
    [t, u] = pde15s(DE, tt, u0, bc, opts);
catch
    errordlg('Error in solution process.', 'chebopbvp error', 'modal');
    varargout{1} = handles;
    return
end

if ( ~guiMode )
    % If we're not running in GUI mode, we can stop here.
    varargout{1} = t;
    varargout{2} = u;
    return
end

%% Post solving steps.

% Store in handles latest solution and time-range (enables exporting later on):
handles.latest.type = 'pde';
handles.latest.solution = u;
handles.latest.solutionT = t;

% Notify the GUI we have a solution available
handles.hasSolution = 1;


if ( guiMode )
    % Switch focus to bottom figure.
    axes(handles.fig_norm)
else
    figure();
end

if ( ~isa(u, 'chebmatrix') )
    waterfall(u, t, 'simple', 'linewidth', defaultLineWidth)
else
    cols = get(0, 'DefaultAxesColorOrder');
    
    % NOTE: Uncomment if we to show legends on the waterfall plot.
%     % Dummy plot to get legends right:
%     for k = 1:size(u,1)
%         plot(0, NaN, 'linewidth', defaultLineWidth, 'color', cols(k, :))
%         hold on
%     end
%     legend(allVarNames);

    % CHEBMATRIX/WATERFALL()
    waterfall(u, t, 'linewidth', defaultLineWidth, 'edgecolors', cols)
end

% Axis labels:
xlabel(indVarName{1})
ylabel(indVarName{2})
zlabel(allVarNames)

% Much pretty. Wow.
view([322.5 30])
box off
grid on    

% Output handles:
varargout{1} = handles;

end