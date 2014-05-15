function varargout = solveguipde(guifile, handles)

% TODO:  Documentation.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/chebfun/ for Chebfun information.

% Create a domain and the linear function on that domain. We use xt for the
% linear function, later in the code we will be able to determine whether x
% or t is used for the linear function.

% Handles will be an empty variable if we are solving without using the GUI
if ( nargin < 2 )
    guiMode = 0;
else
    guiMode = 1;
end

opts = pdeset;
defaultTol = opts.Eps;
defaultLineWidth = 2;

if ( guiMode )
    set(handles.fig_sol, 'Visible', 'On');
    set(handles.fig_norm, 'Visible', 'On');
    cla(handles.fig_sol, 'reset')
    cla(handles.fig_norm, 'reset')
    handles.gui = 1;
else
    handles.gui = 0;
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
% linear function in the problem.
[deString allVarString indVarNameDE pdeVarName pdeflag ignored allVarNames] ...
    = setupFields(guifile, deInput, 'DE');
handles.varnames = allVarNames;
if ( ~any(pdeflag) )
    s = ['Input does not appear to be a PDE, or at least is not a ' ...
         'supported type. Perhaps you need to switch to ''ODE'' mode?'];
    error('Chebgui:SolveGUIpde', s);
end
idx = strfind(deString, ')');

scalarProblem = length(allVarNames) == 1;

% Obtain the independent variable name appearing in the initial condition
if ( ~isempty(initInput{1}) )
    % If we have a scalar problem, we're OK with no dependent variables
    % appearing in the initial guess
    if ( scalarProblem ) 
        [initString ignored indVarNameInit] = ...
            setupFields(guifile, initInput, 'INITSCALAR', allVarString);
        % If the initial guess was just passed a constant, indVarNameInit will
        % be empty. For consistency (allowing us to do the if-statement below),
        % convert to an empty cell.
        if ( isempty(indVarNameInit) )
            indVarNameInit = {''};
        end       
    else
        [initString ignored indVarNameInit] = ...
            setupFields(guifile, initInput, 'INIT', allVarString);        
    end
else
    indVarNameInit = {''};
end

% Make sure we don't have a disrepency in indVarNames. Create a new
% variable indVarName which contains both independent variable, the first
% entry corresponds to space, the second to time.
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
    % TODO:  If never empty for a PDE, why not raise an error?
    indVarName{2} = indVarNameDE{2}; % Should never be empty for a PDE though.
else
    indVarName{2} = 't'; % Default value
end

if ( strcmp(indVarName{1}, indVarName{2}) )
     error('Chebgui:SolveGUIpde', ...
        'The same variable appears to be used as space and time variable');
end
handles.indVarName = indVarName;
opts.handles = handles;

% Support for sum and cumsum
if ( ~isempty(strfind(deString(idx(1):end), 'cumsum(')) )
    sops = {',sum,cumsum'};
elseif ( ~isempty(strfind(deString(idx(1):end), 'sum(')) )
    sops = {',sum'};
else
    sops = {''};
end

% Create a string with the variables used in the problem
variableString = [',', indVarName{2}, ',', indVarName{1}];

deString = [deString(1:idx(1)-1), variableString, sops{:}, ...
    deString(idx(1):end)];
%     deString = strrep(deString,'diff','D');

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

if ( ~isempty(lbcInput{1}) )
    lbcString = setupFields(guifile, lbcInput, 'BC', allVarString);
    idx = strfind(lbcString, ')');
    if ( ~isempty(idx) )
        
        % Support for sum and cumsum
        if ( ~isempty(strfind(lbcString(idx(1):end), 'cumsum(')) )
            sops = {',sum,cumsum'};
        elseif ( ~isempty(strfind(lbcString(idx(1):end),'sum(')) )
            sops = {',sum'};
        else
            sops = {''};
        end
        
        lbcString = [lbcString(1:idx(1)-1), variableString, 'diff', ...
            sops{:}, lbcString(idx(1):end)];
        %             lbcString = strrep(lbcString,'diff','D');
    end
    LBC = eval(lbcString);
else
    LBC = [];
end

if ( ~isempty(rbcInput{1}) )
    rbcString = setupFields(guifile, rbcInput, 'BC', allVarString);
    idx = strfind(rbcString, ')');

    if ( ~isempty(idx) )
        % Support for sum and cumsum
        if ( ~isempty(strfind(rbcString(idx(1):end), 'cumsum(')) )
            sops = {',sum,cumsum'};
        elseif ( ~isempty(strfind(rbcString(idx(1):end), 'sum(')) )
            sops = {',sum'};
        else
            sops = {''};
        end
        
        rbcString = [rbcString(1:idx(1)-1), variableString, 'diff', ...
            sops{:}, rbcString(idx(1):end)];
        %             rbcString = strrep(rbcString,'diff','D');
    end
    RBC = eval(rbcString);
else
    RBC = [];
end

if ( isempty(lbcInput) && isempty(rbcInput) )
    error('chebfun:bvpgui','No boundary conditions specified');
end

if ( isempty(tolInput) )
    tolNum = defaultTol;
else
    tolNum = str2num(tolInput);
end

cpref = chebfunpref();
if ( tolNum < cpref.eps )
    warndlg('Tolerance specified is less than current chebfun epsilon', ...
        'Warning','modal');
    uiwait(gcf)
end

% Boundary condition
if ( periodic )
    bc = 'periodic';
else
    bc = [];
    bc.left = LBC;
    bc.right = RBC;
end

tol = tolNum;

% Set up the initial condition

if ( iscellstr(initInput) )
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
    lenu0 = 0;
    
    if ( isempty(order) && scalarProblem )
        u0 = chebfun(vectorize(initInput{1}), dom);
        u0 = simplify(u0, tol);
        lenu0 = length(u0);
    else
        for k = 1:length(inits)
            initLoc = find(order == k);
            init_k = vectorize(inits{initLoc});
            u0k = chebfun(init_k, dom);
            u0k = simplify(u0k, tol);
            u0 = [u0, u0k];
            lenu0 = max(lenu0, length(u0k));
        end
    end
else
    initInput = vectorize(initInput);
    equalSign = find(initInput == '=');
    if ( isempty(equalSign) )
        equalSign = 0;
    end
    initInput = initInput(equalSign+1:end);
    u0 =  chebfun(initInput, dom);
    u0 = simplify(u0, tol);
    lenu0 = length(u0);
end

% TODO:  Delete if no longer needed.
% if ischar(initInput)
%     initInput = vectorize(initInput);
%     u0 =  chebfun(initInput,[a b]);
%     u0 = simplify(u0,tol);
%     lenu0 = length(u0);
% else
%     u0 = chebfun;
%     lenu0 = 0;
%     for k = 1:numel(initInput)
%         guess_k = vectorize(initInput{k});
%         u0k = chebfun(guess_k,[a b]);
%         u0k = simplify(u0k,tol);
%         u0(:,k) =  u0k;
%         lenu0 = max(lenu0,length(u0k));
%     end
% end

% gather options
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
ylim1 = guifile.options.fixYaxisLower;
ylim2 = guifile.options.fixYaxisUpper;
if ( ~isempty(ylim1) && ~isempty(ylim2) )
    opts.YLim = [str2num(ylim1) str2num(ylim2)];
end
opts.PlotStyle = {'LineWidth', defaultLineWidth};
% TODO:  Delete if no longer needed.
% plotstyle = get(handles.input_plotstyle,'String');
% if ~isempty(plotstyle)
%     opts.PlotStyle = [opts.PlotStyle ',' plotstyle] ;
% end

% Options for fixed N
if ( ~isempty(guifile.options.fixN) )
    opts.N = str2num(guifile.options.fixN);
end

% TODO:  Delete if no longer needed.
% if ~iscell(pdeVarName)
%     pdeVarName = {pdeVarName};
% end
% k = 0;
% idx = [];
% while isempty(idx)
%     k = k+1;
%     idx = strfind(pdeVarName{k},'_');
% end
% timeVarName = pdeVarName{k}((idx+1):end);
% handles.indVarName = {indVarName{1},timeVarName};
% indVarName{2} = timeVarName;

try
    [t u] = pde15s(DE, tt, u0, bc, opts);
catch ME
    errordlg('Error in solution process.', 'chebopbvp error', 'modal');
    return
end

if ( ~guiMode )
    varargout{1} = t;
    varargout{2} = u;
end

% Store in handles latest chebop, solution, vector of norm of updates etc.
% (enables exporting later on)
handles.latest.type = 'pde';
handles.latest.solution = u;
handles.latest.solutionT = t;
% Notify the GUI we have a solution available
handles.hasSolution = 1;

if ( guiMode )
    axes(handles.fig_norm)
else
    figure();
end

if ( ~iscell(u) )
    %     surf(u,t,'facecolor','interp')
    waterfall(u, t, 'simple', 'linewidth', defaultLineWidth)
    xlabel(indVarName{1}), ylabel(indVarName{2}), zlabel(allVarNames)
else
    cols = get(0, 'DefaultAxesColorOrder');

    for k = 1:numel(u)
        plot(0, NaN, 'linewidth', defaultLineWidth, 'color', cols(k, :))
        hold on
    end

    legend(allVarNames);

    for k = 1:numel(u)
        waterfall(u{k}, t, 'simple', 'linewidth', defaultLineWidth, ...
            'edgecolor', cols(k, :))
        hold on
        xlabel(indVarName{1})
        ylabel(indVarName{2})
    end

    view([322.5 30])
    box off
    grid on
    
    hold off
end

varargout{1} = handles;
