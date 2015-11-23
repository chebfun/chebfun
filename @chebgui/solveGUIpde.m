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

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Handles will be an empty variable if we are solving without using the GUI
if ( nargin < 2 )
    guiMode = 0;
else
    guiMode = 1;
end

% Call the exportInfo method of the chebguiExporterEIG class, which takes care
% of extracting most information we need from GUIFILE.
expInfo = chebguiExporterPDE.exportInfo(guifile);

% Extract information from the EXPINFO struct.
dom = str2num(expInfo.dom);
tt = str2num(expInfo.tt);
deString = expInfo.deString;
allVarNames = expInfo.allVarNames;
initInput = expInfo.initInput;
indVarName = expInfo.indVarName;
pdeflag = expInfo.pdeflag;
periodic = expInfo.periodic;
lbcString = expInfo.lbcString;
rbcString = expInfo.rbcString;
pdeSolver = eval(['@', expInfo.pdeSolver]);

% Check that we don't have any breakpoints.
if ( length(dom) > 2 )
    warning('CHEBFUN:CHEBGUI:solveGUIpde', ...
        ['PDE solver does not accept domains with breakpoints. ' ...
         'Breakpoints are being ignored.']);
    dom = dom([1 end]);
end


% Are we actually trying to solve a PDE?
if ( ~any(pdeflag) )
    s = ['Input does not appear to be a PDE, or at least is not a ' ...
         'supported type. Perhaps you need to switch to ''ODE'' mode?'];
    error('CHEBFUN:CHEBGUI:solveGUIpde:solveGUIpde', s);
end

% Are we solving a scalar problem?
scalarProblem = length(allVarNames) == 1;

% Store useful variable names in HANDLES:
handles.indVarName = indVarName;
handles.varnames = allVarNames;

% Convert the DE string to a proper anonymous function using eval:
DE = eval(deString);

% Do we have a left BC specified?
if ( ~isempty(lbcString) )
    LBC = eval(lbcString);
else
    LBC = [];
end

% Do we have a right BC specified?
if ( ~isempty(rbcString) )    
    RBC = eval(rbcString);
else
    RBC = [];
end


% Setup up boundary condition to be passed to the pde15s() method:
if ( periodic )
    bc = 'periodic';
elseif ( isempty(lbcString) && isempty(rbcString) )
    % No boundary conditions are a no-go:
    error('CHEBFUN:CHEBGUI:solveGUIpde:bvpgui', ...
        'No boundary conditions specified');
else
    bc.left = LBC;
    bc.right = RBC;
end

% Start gathering options. Start by obtaining a PDE options struct.
opts = pdeset;
defaultTol = opts.Eps;
defaultLineWidth = 2;

% Check that we're not imposing too strict tolerance:
tolInput = guifile.tol;
if ( isempty(tolInput) )
    tol = defaultTol;
else
    tol = str2double(tolInput);
end

% Find out what the current CHEBFUN eps is set at.
cpref = chebfunpref();
if ( tol < cpref.eps )
    warndlg('Tolerance specified is less than current chebfun epsilon', ...
        'Warning','modal');
    uiwait(gcf)
end

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
    else
        for k = 1:length(inits)
            initLoc = find(order == k);
            init_k = vectorize(inits{initLoc});
            u0k = chebfun(init_k, dom);
            u0 = [u0, u0k];
        end
    end
else
    % Only one line specifying the initial condition (for a scalar problem):
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
    set(handles.fig_sol, 'fontsize', handles.fontsizePanels);
    set(handles.fig_norm, 'fontsize', handles.fontsizePanels);
    handles.gui = 1;
else
    handles.gui = 0;
end

% Continue setting up options.

% Give the OPTS structure the handle to the chebguiwindow figure, so that pde15s
% can plot to it:
opts.handles = handles;

opts.Eps = tol;
if ( ~all(pdeflag) )
    opts.PDEflag = pdeflag;
end

% Plot during PDE solving?
opts.Plot = expInfo.doplot;
opts.HoldPlot = expInfo.dohold;

% Do we want to fix the y-limits while solving the PDE?
if ( ~isempty(expInfo.ylim1) && ~isempty(expInfo.ylim2) )
    opts.YLim = [str2double(expInfo.ylim1), str2double(expInfo.ylim2)];
end

% We want thicker lines on the plot:
opts.PlotStyle = {'LineWidth', defaultLineWidth};

% Options for fixed N:
if ( ~isempty(expInfo.fixN) )
    opts.N = str2double(expInfo.fixN);
end

% Try solving the PDE!
try
    [t, u] = pdeSolver(DE, tt, u0, bc, opts);
catch ME
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
    waterfall(u, t)
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
    waterfall(u, t,'edgecolors', cols)
end

% Update the fontsize of the bottom plot 
set(handles.fig_norm, 'fontsize', handles.fontsizePanels);
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
