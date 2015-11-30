function varargout = pdeSolve(pdeFun, tIn, u0, bc, varargin)
%PDESOLVE   Solve PDEs using Chebfun.
%
%   PDESOLVE() solves solves an initial-boundary value problem via a method of
%   lines approach. See @CHEBFUN/PDE15S() for documentation.
%
% See also PDE15S, PDE23T, PDESET.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

global DIFFORDER SYSSIZE
DIFFORDER = 0;
SYSSIZE = 0;

% Default options:
tol = 1e-6;             % 'eps' in Chebfun terminology
% The default behaviour with no outputs is to plot. If the method is called with
% outputs, by default, we don't plot.
doPlot = ( nargout == 0 );
doHold = false;         % Hold plot?
plotOpts = {'-'};       % Plotting style
adjustBCs = true;       % Adjust inconsistent BCs
throwBCwarning = true;  % Throw a warning for inconsistent BCs
timeChunks = 51;        % Default number of time slices if not specified

% Parse the variable inputs:
if ( numel(varargin) == 2 )
    opt = varargin{1};
    opt.N = varargin{2};
elseif ( numel(varargin) == 1 )
    if ( isstruct(varargin{1}) )
        opt = varargin{1};
    else
        opt = pdeset;
        opt.N = varargin{1};
    end
else
    opt = pdeset;
end
optN = opt.N;
if ( isempty(optN) )
    optN = NaN;
end

% Get the domain:
DOMAIN = domain(u0, 'ends');

if ( isfield(opt, 'AdjustBCs') && ~isempty(opt.AdjustBCs) && ~opt.AdjustBCs )
    throwBCwarning = false;
    adjustBCs = false;
end

% PDE solver options:
if ( ~isempty(opt.Eps) )
    tol = opt.Eps;
end
if ( ~isempty(opt.Plot) )
    doPlot = strcmpi(opt.Plot, 'on');
end
if ( ~isempty(opt.HoldPlot) )
    doHold = strcmpi(opt.HoldPlot, 'on') || strcmp(opt.HoldPlot, '1');
end
if ( ~isempty(opt.PlotStyle) )
    plotOpts = opt.PlotStyle;
end
if ( ~isempty(opt.ODESolver) )
    ODESOLVER = opt.ODESolver;
else
    ODESOLVER = @ode15s;
end

% Set time chunks:
if ( length(tIn) == 2 )
    tt = linspace(tIn(1), tIn(2), timeChunks);
else
    tt = tIn;
end

userMassSet = false;
ISPERIODIC = false;
DONE = false;
COUNTER = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%  PARSE INPUTS TO PDEFUN  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the size of the system, i.e., number of dependent variables.
SYSSIZE = min(size(u0));
[pdeFun, varNamesParsed] = parseFun(pdeFun);
if ( isfield(opt, 'difforder') )
    DIFFORDER = opt.difforder;
else
    getDIFFORDER(pdeFun);
end

if ( (max(DIFFORDER) < 2) && isequal(ODESOLVER, @ode15s) )
    warning('CHEBFUN:PdeSolver', ...
        'PDE23T() is recommended for non-diffusive problems.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  EVENT SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ~isnan(optN) )
    opt.OutputFcn = @nonAdaptiveEvent;
else
    opt.OutputFcn = @adaptiveEvent;
end

    function status = nonAdaptiveEvent(t, U, flag)
        % This event is called at the end of each chunk in non-adaptive mode.
        status = false;
        if ( ~isempty(flag) )
            return
        end
        
        % Sometimes we get given more than one time slice.
        for kk = 1:numel(t)
            
            if ( ~any(abs(tSpan - t(kk)) < 1e-6) )
                % This is not a designated time slice!
                continue
            end
                
            % Reshape solution:
            Uk = reshape(U(:,kk), n, SYSSIZE);
            uCurrent = chebfun(Uk, DOMAIN, 'tech', techHandle);
            tCurrent = t(kk);
            % Store for output:
            COUNTER = COUNTER + 1;
            uOut{COUNTER} = uCurrent;

            % Plot current solution:
            if ( doPlot )
                plotFun(uCurrent, t(kk));
            end
        end
        
        if ( guiFlag )
            status = guiEvent(status);
        end

    end

    function status = adaptiveEvent(t, U, flag)
        % This event is called at the end of each chunk in adaptive mode.

        status = false;
        if ( ~isempty(flag) )
            return
        end

        % Sometimes we get given more than one time slice.
        for kk = 1:numel(t)
            % Reshape solution:
            Uk = reshape(U(:,kk), currentLength, SYSSIZE);

            % Happiness check:
            c = (1+sin(1:SYSSIZE)).'; % Arbitrarily linear combination.
            Uk2 = (Uk*c/sum(c));
            uk2 = tech.make(Uk2, pref);
            [ishappy, cutoff] = happinessCheck(uk2, Uk2, [], [], pref);

            if ( ishappy )  

                if ( ~any(abs(tSpan - t(kk)) < 1e-6) )
                    % This is not a designated time slice!
                    continue
                end

                % Store these values:
                tCurrent = t(kk);
                uCurrent = chebfun(Uk, DOMAIN, 'tech', techHandle);
                uCurrent = simplify(uCurrent);
                
                COUNTER = COUNTER + 1;
                uOut{COUNTER} = uCurrent;
                tOut(COUNTER) = tCurrent;

                % Plot current solution (if plotting is on)
                if ( doPlot )
                    plotFun(uCurrent, tCurrent);
                end
                % TODO: Re-insert this?
%                 % If we have 2.5 times as many coefficients as we need then 
%                 % shorten the representation and cause the integrator to stop. 
%                 if ( cutoff < 0.4*n && n > 17)
%                     currentLength = round(1.25*cutoff)';
%                     %currentLength = floor( currentLength / 1.5  );
%                     currentLength = currentLength + 1 - rem(currentLength,2);
%                     currentLength = max(currentLength, 17);
%                     status = true;
%                     return
%                 end
                
            else 

                % Increase length and bail out:
                currentLength = 2*currentLength-1;
                status = true;
                break
                
            end        
            
            if ( guiFlag )
                status = guiEvent(status);
            end
            
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOTTING SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine which figure to plot to (for CHEBGUI) and set default display values
% for variables.
YLim = opt.YLim;
gridOn = 0;
guiFlag = false;
if ( isfield(opt, 'handles') )
    if ( opt.handles.gui )
        guiFlag = true;
        axesSol = opt.handles.fig_sol;
        axesNorm = opt.handles.fig_norm;
        gridOn = opt.handles.guifile.options.grid;
        solveButton = opt.handles.button_solve;
        clearButton = opt.handles.button_clear;
        panelSol = opt.handles.panel_figSol;
    end
    varNames = opt.handles.varnames;
    xLabel = opt.handles.indVarName{1};
    tlabel = opt.handles.indVarName{2};
    fontsize = opt.handles.fontsizePanels;
    set(axesSol, 'fontsize', fontsize);
else
    % Obtain the correct variable names from the parsing of the operator
    if ( length(varNamesParsed) == numColumns(u0) )
        % Space and time variables did not get passed explicitly to operator
        tlabel = 't';
        xLabel = 'x';
    else
        tlabel = varNamesParsed{1};
        xLabel = varNamesParsed{2};
    end
    varNames = varNamesParsed(end-numColumns(u0)+1:end);
    if ( doPlot )    % Set up ploting in non-gui mode
        axesSol = gca;
        cla(axesSol)
        xlabel(axesSol, xLabel);
        % Bring attention to the figure
        shg
    end
end

% Set plot options
if ( doPlot )
    set(axesSol, 'NextPlot','replacechildren');
    
    % Fix x limits
    set(axesSol, 'xLim', u0.domain);
    
    % Fix y limits if they were specified
    if ( ~isempty(YLim) )
        set(axesSol, 'ylim', YLim);
    end
    
    % Initialize lines for plotting
    hLines = plot(axesSol, DOMAIN, NaN(length(DOMAIN), numColumns(u0)), ...
        plotOpts{:});
    
    % Grid on?
    if ( gridOn )
        grid(axesSol, 'on')
    else
        grid(axesSol, 'off')
    end
    
    % Hold on?
    if ( doHold )
        hold(axesSol, 'on')
    else
        hold(axesSol, 'off')
    end
    
    % For plotting, it's useful to know whether we're running in old or new
    % Matlab graphics mode
    if ( ~verLessThan('matlab', '8.4') )
        newMatlabVersion = true;
    else
        newMatlabVersion = false;
    end
end


    function status = guiEvent(status)
        %GUIEVENT   Deal with GUI events ('stop', 'pause', etc). 
        % OUTPUTS:
        %   status = true exits the current time chunk.
        %   DONE = true exits PDE solver. (Note DONE is a GLOBAL variable).

        % Interrupt computation if stop or pause button is pressed in the GUI.
        if ( strcmp(get(solveButton, 'String'), 'Solve') )
            % Stop.
            tt = tt( tt <= tCurrent );
            uOut(COUNTER+1:end) = [];
            status = true;
            DONE = true;
        elseif ( strcmp(get(clearButton, 'String'), 'Continue') )
            % Wait, pause.
            
            % Plot a waterfall plot to the bottom figure window:
            axes(axesNorm)
            uuTmp = prepare4output(uOut(1:COUNTER));
            if ( SYSSIZE == 1 )
                waterfall(uuTmp, tt(tt<=tCurrent), 'linewidth', 2);
            else
                cols = get(0, 'DefaultAxesColorOrder');
                waterfall(uuTmp, tt(tt<=tCurrent), 'linewidth', 2, ...
                    'EdgeColors', cols);
            end
            xlabel(xLabel), ylabel(tlabel), zlabel(varNames)
            if ( gridOn )
                grid on
            end
            view([322.5 30])
            box off
            axes(axesSol)
            % Hang around until 'CONTINUE' or 'STOP' is presed.
            waitfor(clearButton, 'String');
            % Call again to see if 'STOP' was pressed.
            status = guiEvent(status);
            
            set(axesNorm, 'fontsize', fontsize);
        end
    end

    function plotFun(U, t)
        %PLOTFUN    Plot current solution U at a time t.

        % Obtain plotting data for the solution
        Udata = plotData(U);
        UdataX = Udata.xLine;
        UdataY = Udata.yLine;
        
        % Update the YData of the lines plotted
        if ( ~doHold )
            for hCounter = 1:length(hLines)
                set(hLines(hCounter), 'XData', UdataX);
                set(hLines(hCounter), 'YData', UdataY(:, hCounter)');
            end
        else
            % Reset color cycle prior to point plot if running R2014b.
            if ( newMatlabVersion )
                set(axesSol, 'ColorOrderIndex', 1);
            end
            hLines = plot(axesSol, UdataX, UdataY, plotOpts{:});
        end
        
        % Update the title, either of the GUI panel or the figure
        if ( guiFlag )
            set(panelSol, 'Title', sprintf('%s = %.3f', tlabel, t))
        elseif ( nargin > 1 )
            title(sprintf('%s = %.3f', tlabel, t))
        end
        
        % Update the plot
        drawnow

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%% PARSE BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ischar(bc) && (strcmpi(bc, 'neumann') || strcmpi(bc, 'dirichlet')) )
    bc = struct( 'left', bc, 'right', bc);
elseif ( iscell(bc) && numel(bc) == 2 )
    bc = struct( 'left', bc{1}, 'right', bc{2});
end
if ( isstruct(bc) && ~isfield(bc, 'middle') )
    bc.middle.op = [];
end

% Initialise some global variables.
leftNonlinBCLocs = [];
middleNonlinBCLocs = [];
rightNonlinBCLocs = [];
BCRHS = {};

if ( ischar(bc) && any(strcmpi(bc, {'periodic', 'trig'})) )
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% PERIODIC BCS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ISPERIODIC = true;
    
    % One can still use a Chebyshev basis and enforce periodic conditions by
    % enforcing suitable constraints on the derivative. However, using a
    % periodic basis (for now, trigtech) is much more efficient.
%     r = cell(sum(DIFFORDER), 1);
%     count = 1;
%     for j = 1:SYSSIZE
%         for k = 0:DIFFORDER(j)-1
%             c = (diff(DOMAIN)/2)^k;
%             A = @(n) [1 zeros(1, n-2) -1]*chebcolloc2.diffmat(n, k)*c;
%             r{count} = @(n) [zeros(1, (j-1)*n) A(n) zeros(1,(SYSSIZE-j)*n)];
%             count = count + 1;
%         end
%     end
%     bc = struct( 'left', [], 'right', []);
%     bc.left.op = r(1:2:end);
%     bc.right.op = r(2:2:end);
%     BCRHS = num2cell(zeros(1, numel(r)));
    
else
    %% %%%%%%%%%%%%%%%%%%%%%%%%% NONPERIODIC BCS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ( isfield(bc, 'left') && ~isfield(bc, 'right') )
        bc.right = [];
    elseif ( isfield(bc, 'right') && ~isfield(bc, 'left') )
        bc.left = [];
    elseif ( ~isfield(bc, 'left') && ~isfield(bc, 'right') )
        bc = struct('left', [], 'right', [], 'middle', bc);
    end
    
    % Deal with struct and numeric input:
    bc.left = dealWithStructInput(bc.left);
    bc.right = dealWithStructInput(bc.right);
    
    if ( isempty(bc.left) )
        bc.left = struct('op', []);
    elseif ( ischar(bc.left) || (iscell(bc.left) && ischar(bc.left{1})) )
        %% %%%%%%%%%%%%%%%%%%%%% DIRICHLET AND NEUMANN BCS (LEFT) %%%%%%%%%%%%%%
        if ( iscell(bc.left) )
            v = bc.left{2};
            bc.left = bc.left{1};
        else
            v = 0;
        end
        if ( ~isnumeric(v) )
            error('CHEBFUN:CHEBFUN:pde15s:nonNumericVal1', ...
                'For BCs of the form {char, val} val must be numeric.')
        end
        if ( strcmpi(bc.left, 'dirichlet') )
            A = @(n) [1 zeros(1, n - 1)];
        elseif ( strcmpi(bc.left, 'neumann') )
            % TODO: Make left diff operator explicitly.
            A = @(n) [1 zeros(1, n-1)]*chebcolloc2.diffmat(n)*diff(DOMAIN)/2;
        else
            error('CHEBFUN:CHEBFUN:pde15s:bcSyntax1', 'Unknown BC syntax');
        end
        bc.left = struct('op', []);
        bc.left.op = cell(SYSSIZE, 1);
        for k = 1:SYSSIZE
            bc.left.op{k} = @(n) [zeros(1, ( k -1)*n) A(n) ...
                zeros(1 , (SYSSIZE - k)*n)];
        end
        BCRHS = num2cell(repmat(v, SYSSIZE, 1));
    elseif ( numel(bc.left) == 1 && isa(bc.left, 'function_handle') )
        %% %%%%%%%%%%%%%%%%%%%%%%%%%  GENERAL BCS (LEFT)  %%%%%%%%%%%%%%%%%%%%%%
        op = parseFun(bc.left, 'lbc');
        uTmp = chebdouble(ones(1, SYSSIZE));
        sizeOp = size(op(0, mean(DOMAIN), uTmp));
        leftNonlinBCLocs = 1:max(sizeOp);
        bc.left = struct('op', []);
        bc.left.op = {@(n) zeros(max(sizeOp), SYSSIZE*n)}; % Dummy entries.
        BCRHS = num2cell(zeros(1, max(sizeOp)));
        leftNonlinBCFuns = op;
    else
        error('CHEBFUN:CHEBFUN:pde15s:bcSyntax2', 'Unknown BC syntax');
    end
    
    if ( isfield(bc, 'middle') && isa(bc.middle, 'function_handle') )
        %% %%%%%%%%%%%%%%%%%%%%%  GENERAL BCS (MIDDLE)  %%%%%%%%%%%%%%%%%%%%%%%%
        op = parseFun(bc.middle, 'bc');
        uTmp = chebdouble(ones(1, SYSSIZE));
        sizeOp = size(op(0, mean(DOMAIN), uTmp));
        middleNonlinBCLocs = 1:max(sizeOp);
        bc.middle = [];
        bc.middle.op = {@(n) zeros(max(sizeOp), SYSSIZE*n)}; % Dummy entries.
        BCRHS = num2cell(zeros(1, max(sizeOp)));
        middleNonlinBCFuns = op;
    end
    
    if ( isempty(bc.right) )
        bc.right = struct('op', []);
    elseif ( ischar(bc.right) || (iscell(bc.right) && ischar(bc.right{1})) )
        %% %%%%%%%%%%%%%%%%%%%%% DIRICHLET AND NEUMANN BCS (RIGHT) %%%%%%%%%%%%%
        if ( iscell(bc.right) )
            v = bc.right{2};
            bc.right = bc.right{1};
        else
            v = 0;
        end
        if ( ~isnumeric(v) )
            error('CHEBFUN:CHEBFUN:pde15s:nonNumericVal1', ...
                'For BCs of the form {char, val} val must be numeric.')
        end
        if ( strcmpi(bc.right, 'dirichlet') )
            A = @(n) [zeros(1, n-1), 1];
        elseif ( strcmpi(bc.right, 'neumann') )
            % TODO: Make right diff operator explicitly.
            A = @(n) [zeros(1, n-1) 1]*chebcolloc2.diffmat(n)*diff(DOMAIN)/2;
        else
            error('CHEBFUN:CHEBFUN:pde15s:bcSyntax3', 'Unknown BC syntax');
        end
        bc.right = struct('op', []);
        bc.right.op = cell(SYSSIZE, 1);
        for k = 1:SYSSIZE
            bc.right.op{k} = @(n) [zeros(1,(k-1)*n) A(n) zeros(1,(SYSSIZE-k)*n)];
        end
        BCRHS = [BCRHS num2cell(repmat(v, SYSSIZE, 1))];
        
    elseif ( numel(bc.right) == 1 && isa(bc.right, 'function_handle') )
        %% %%%%%%%%%%%%%%%%%%%%%%%%  GENERAL BCS (RIGHT)  %%%%%%%%%%%%%%%%%%%%%%
        op = parseFun(bc.right, 'rbc');
        uTmp = chebdouble(ones(1, SYSSIZE));
        sizeOp = size(op(0, mean(DOMAIN), uTmp));
        rightNonlinBCLocs = 1:max(sizeOp);
        bc.right = struct('op', []);
        bc.right.op = {@(n) zeros(max(sizeOp), SYSSIZE*n)};
        BCRHS = [BCRHS num2cell(zeros(1, max(sizeOp)))];
        rightNonlinBCFuns = op;
    else
        error('CHEBFUN:CHEBFUN:pde15s:bcSyntax4', 'Unknown BC syntax');
    end
    
end

if ( isstruct(bc) && ~isfield(bc, 'middle') )
    bc.middle.op = [];
end

%% %%%%%%%%%%%%%%%%%%%%%%% SUPPORT FOR COUPLED BVP-PDES! %%%%%%%%%%%%%%%%%%%%%%%
% Experimental feature for coupled ode/pde systems: (An entry equal to 1 denotes
% that the corresponding variable appears with a time derivative. 0 otherwise.)
if ( isfield(opt, 'PDEflag') && ~isempty(opt.PDEflag) )
    pdeFlag = opt.PDEflag;
else
    pdeFlag = true;
end
if ( numel(pdeFlag) == 1 )
    pdeFlag = repmat(pdeFlag, 1, SYSSIZE);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERIODIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ~ISPERIODIC )
    techHandle = @chebtech2;
    points = @chebpts;
    mydouble = @chebdouble;
else
    techHandle = @trigtech;
    points = @trigpts;
    mydouble = @trigdouble;  
    u0 = chebfun(u0, 'trig', 'eps', tol);
    u0 = simplify(u0, tol);
end
tech = techHandle();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% SETUP TOLERANCES AND INITIAL CONDITION %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ODE tolerances:
aTol = odeget(opt, 'AbsTol', tol);
rTol = odeget(opt, 'RelTol', tol);
if ( isnan(optN) )
    aTol = min(aTol, tol);
    rTol = min(rTol, tol);
end
opt.AbsTol = aTol;
opt.RelTol = rTol;

% Check for (and try to remove) piecewise initial conditions:
u0Trans = u0(1).isTransposed;
if ( u0Trans )
    u0 = transpose(u0);
end
for k = 1:numel(u0)
    if ( numel(u0(k).funs) > 1 )
        u0(k) = merge(u0(k), 'all', 1025, tol);
        if ( u0(k).nfuns > 1 )
            error('CHEBFUN:CHEBFUN:pde15s:piecewise', ...
                'Piecewise initial conditions are not supported.');
        end
    end
end

% Simplify initial condition to tolerance or fixed size in optN:
if ( isnan(optN) )
    u0 = simplify(u0);
else
    for k = 1:numel(u0)
        u0(k).funs{1}.onefun = prolong(u0(k).funs{1}.onefun, optN);
    end
end

% Initial condition:
uCurrent = u0;
% Storage:
uOut = cell(1, numel(tt));
uOut{1} = uCurrent;
tOut(1) = tt(1);

% Initialise variables for ONESTEP():
B = []; q = []; rows = []; M = []; P = []; n = [];

% Set the preferences:
pref = tech.techPref();
pref.eps = tol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% TIME CHUNKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot initial condition:
currentLength = length(u0);
if ( doPlot )
    plotFun(u0, tt(1));
    leg = legend(axesSol, varNames);
    set(leg, 'HandleVisibility', 'off')
    drawnow
end
DONE = false;
if ( ~isnan(optN) )
    % Non-adaptive in space:
    tSpan = tt;
    x = points(optN, DOMAIN);
    solvePDE(tSpan); % Do all chunks at once!
else
    % Adaptive in space
    
    % Min-length is 9
    currentLength = max(currentLength, 9);
    tCurrent = tt(1);
    tSpan = tt;
    while ( tCurrent < tt(end) && ~DONE )
        tSpan(tSpan < tCurrent) = [];
        x = points(currentLength, DOMAIN);
        solvePDE(tSpan);
    end

end

if ( doPlot && ~doHold )
    hold off
end

% Ensure tOut is a column vector
tOut = tOut(:);

uOut = prepare4output(uOut);

switch nargout
    case 0
    case 1
        if ( length(tIn) == 2 )     % Only T_start and T_end were given.
            uOut = uOut(:,[1,end]); % Strip intermediate steps.
        end
        varargout{1} = uOut;
    case 2
        varargout{1} = tOut;
        varargout{2} = uOut;
    case 1 + size(uOut, 1)
        varargout{1} = tOut;
        for varCounter = 1:size(uOut, 1)
            varargout{varCounter + 1} = chebfun(uOut(varCounter, :));
        end
    otherwise
        error('CHEBFUN:CHEBFUN:pdeSolve:nargout', ...
            'Incorrect number of output arguments.');
end

clear global DIFFORDER
clear global SYSSIZE

    function uOut = prepare4output(uIn)
        % If we only had one dependent variable, return an array valued CHEBFUN
        % instead of a QUASIMATRIX.
        if ( SYSSIZE == 1 )
            uOut = horzcat(uIn{:});
        else
            blocks = cell(SYSSIZE, numel(uIn));
            for kk = 1:SYSSIZE
                blocks(kk,:) = cellfun(@(u) extractColumns(u, kk), uIn, ...
                    'UniformOutput', false);
            end
            uOut = chebmatrix(blocks); % CHEBMATRIX
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ONESTEP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function solvePDE(tSpan)
        % Constructs the result of one time chunk at fixed discretization.

        % Evaluate the chebfun at discrete points:
        U0 = feval(uCurrent, x);

        if ( ISPERIODIC )
            
            % The discretisation length
            n = length(x);
            
        elseif ( isempty(n) || (n ~= length(x)) )
            % This depends only on the size of n. If this is the same, reuse!
            
            % The new discretisation length
            n = length(x);
            
            % Linear constraints:
            bcop = [bc.left.op ; bc.middle.op ; bc.right.op];
            B = cell2mat(cellfun(@(f) feval(f, n), bcop, 'UniformOutput', false));
            
            % Project / mass matrix.
            P = cell(1, SYSSIZE);
            M = cell(1, SYSSIZE);
            
            for kk = 1:SYSSIZE
                xk = chebpts(n-DIFFORDER(kk), DOMAIN, 1);
                P{kk} = barymat(xk, x);
                M{kk} = pdeFlag(kk)*P{kk};
            end
            P = [ 0*B ; blkdiag(P{:})];
            M = [ 0*B ; blkdiag(M{:})];
            rows = 1:size(B, 1);
            
            % Multiply by user-defined mass matrix
            if ( userMassSet )
                M = feval(userMass, n)*M;
            end
            
            % ODE options: (mass matrix)
            opt = odeset(opt, 'Mass', M, 'MassSingular', 'yes', ...
                'MStateDependence', 'none');
            
        end
        
        if ( ~ISPERIODIC )
            % We have to ensure the starting condition satisfies the boundary
            % conditions, or else we will get a singularity in time. We do this by
            % tweaking the BC values to match the reality. Changing the IC itself is
            % much trickier. Find out what the BC deviance from nominal really is:
            BCVALOFFSET = 0;            % recover nominal value in next call
            F = odeFun(tSpan(1),U0(:)); % also assigns to "rows" and "q"

            % If this is for the initial chunk, check whether the initial
            % condition nearly satisfies the BCs.
            % We're quite lax about this, because discretization at low N can
            % cause derivatives to look fairly bad. 
            if ( throwBCwarning && (length(uOut) > 1) && (norm(F(rows)) > 0.05*norm(F)) )
                warning('CHEBFUN:CHEBFUN:pde15s:BadIC',...
                    'Initial state may not satisfy the boundary conditions.')
                throwBCwarning = false;
            end
            if ( adjustBCs )
                BCVALOFFSET = F(rows) - q;
            else
                BCVALOFFSET = 0;
            end
        end
        
        % Solve ODE over time chunk with the selected solver:
        try
            [ignored1, ignored2] = ODESOLVER(@odeFun, tSpan, U0, opt);
        catch ME
            if ( strcmp(ME.identifier, 'MATLAB:odearguments:SizeIC') )
                error('CHEBFUN:CHEBFUN:pde15s:dims', ...
                    'Dimension mismatch. Check boundary conditions.');
            else
                rethrow(ME)
            end
        end 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%  ODEFUN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function F = odeFun(t, U)
            % This is what the ODESOLVER() calls.
            
            % Reshape to n by SYSSIZE:
            U = reshape(U, n, SYSSIZE);
            
            % Evaluate the PDEFUN:
            Utmp = mydouble(U, DOMAIN);
            F = pdeFun(t, x, Utmp);
            F = double(F);
            F = F(:);
            
            if ( ISPERIODIC )
                return
            end
            
            % Enforce boundary constraints:
            
            F = P*F; % Project.
            
            % Get the algebraic right-hand sides: (may be time-dependent)
            for l = 1:numel(BCRHS)
                if ( isa(BCRHS{l}, 'function_handle') )
                    q(l, 1) = feval(BCRHS{l}, t);
                else
                    q(l, 1) = BCRHS{l};
                end
            end
            
            % Replacements for the BC algebraic conditions:
            F(rows) = B*U(:) - q;
            
            % Replacements for the nonlinear BC conditions:
            if ( ~isempty(leftNonlinBCLocs) )
                indx = 1:length(leftNonlinBCLocs);
                uTmp = leftNonlinBCFuns(t, x, Utmp);
                uTmp = double(uTmp);
                if ( size(uTmp, 1) ~= n )
                    uTmp = reshape(uTmp, n, numel(uTmp)/n);
                end
                F(rows(indx)) = uTmp(1, :);
            end
            if ( ~isempty(middleNonlinBCLocs) )
                % TODO: This won't work if there are also left nonlin BCs
                indx = 1:length(middleNonlinBCLocs);
                F(rows(indx)) = double(middleNonlinBCFuns(t, x, Utmp));
            end            
            if ( ~isempty(rightNonlinBCLocs) )
                indx = numel(BCRHS) + 1 - rightNonlinBCLocs;
                uTmp = rightNonlinBCFuns(t, x, Utmp);
                uTmp = double(uTmp);
                if ( size(uTmp, 1) ~= n )
                    uTmp = reshape(uTmp, n, numel(uTmp)/n);
                end
                F(rows(indx)) = fliplr(uTmp(end, :));
            end
            
            % Adjust BC rows by the needed offset from original time:
            F(rows) = F(rows) - BCVALOFFSET;
            
            % Reshape to back to a single column:
            F = F(:);
            
        end
    end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%  MISC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outFun, varNamesParsed] = parseFun(inFun, bcFlag)
% Rewrites the input function handle to call the right DIFF, SUM, methods, etc,
% and convert the input @(t, x, u, v, w, ...) to @(t, x, U).
%
% For PARSEFUN(INFUN), if the INFUN has only one independent variable as an
% input, we assume this is x. PARSEFUN(INFUN, 'lbc') or PARSEFUN(INFUN, 'rbc')
% will assume the variable is t. PARSEFUN(INFUN, 'bc') will assume space, but
% throw a warning.

global SYSSIZE

% Ensure backwards compatibility by processing the input function.
[inFun, varNamesParsed] = backCompat(inFun);

% V4: We don't accept only time or space as input args (i.e., both or nothing).
% V5: Actually, we now accept @(u, x) diff(u, 2) + sin(x).*u (for CHEBOPs).
Nind = nargin(inFun) - SYSSIZE;
if ( Nind == 0 )
    outFun = @(t, x, u) conv2cell(inFun, u);
elseif ( Nind == 1 )
    if ( nargin > 1 )
        if ( strcmp(bcFlag, 'bc') )
            warning('CHEBFUN:CHEBFUN:pde15s:bcInput', ...
                ['Please input both time and space independent variable to ' ...
                 '.BC function handle inputs.\nAssuming given variable is ' ...
                 'time.']);
             outFun = @(t, x, u) conv2cell(inFun, u, x);
        else
            outFun = @(t, x, u) conv2cell(inFun, u, t);
        end
    else
        outFun = @(t, x, u) conv2cell(inFun, u, x);
    end
elseif ( Nind == 2 )
    outFun = @(t, x, u) conv2cell(inFun, u, t, x);
else
    error('CHEBFUN:CHEBFUN:pde15s:inputs_ind', ...
        ['Incorrect number of independent variables in input function. ' ...
         '(Must be 0, 1, or 2).']);
end

end

function newFun = conv2cell(oldFun, u, varargin)
% This function allows the use of different variables in the anonymous function,
% rather than using the clunky quasi-matrix notation.
%
% Inputs: 
%   oldFun in the original function handle
%   u will be the solution we evaluate at
%   varargin possibly contains the time and space variables
% Output
%   newFun the modified function handle

global SYSSIZE
% Note. This is faster than mat2cell!
tmpCell = cell(1, SYSSIZE);
for qk = 1:SYSSIZE
    tmpCell{qk} = extractColumns(u, qk);
end
newFun = oldFun(varargin{:}, tmpCell{:});
end

function getDIFFORDER(pdeFun)
% Extract the DIFFORDER by evaluating the operator at some NaN values.
global DIFFORDER SYSSIZE
% Find the order of each of the variables:
tmp = chebdouble(zeros(SYSSIZE));
v = pdeFun(0, 0, tmp);
DIFFORDER = get(v, 'diffOrder');
end

function out = dealWithStructInput(in)
% Parse the inputs for BCs.
% Throw a warning or an error for a struct. Convert numeric value to 'dirichet'.
if ( isstruct(in) )
    if ( numel(in) == 1)
        warning('CHEBFUN:CHEBFUN:pde15s:bcstruct', ...
            'PDE15S no longer supports struct inputs for bc.left and bc.right.')
        if ( isfield(in, 'val') )
            out = {in.op, in.val};
        else
            out = in.op;
        end
    else
        error('CHEBFUN:CHEBFUN:pde15s:bcstruct', ...
            'PDE15S no longer supports struct inputs for bc.left and bc.right.')
    end
elseif ( isempty(in) )
    out = [];
elseif ( isnumeric(in) )
    out = {'dirichlet', in};
else
    out = in;
end
end

function [outFun, varNamesOut] = backCompat(inFun)
% In V4 PDE15S required the user to pass in dummy function handles for the
% differential operators. For example, u'' + u' + sum(u) would have needed
%  pdefun = @(u, t, x, diff, sum) diff(u, 2) + diff(u) + sum(u).
% V5 deals with this in a different way (by using the CHEBDOUBLE class), but we
% would still like to support the old syntax (at least for now).
%
% To do this, we parse the input function for the strings 'diff', 'Diff', and
% 'D'. If we find any of these, we assume the old syntax is being used, and
% redefine the function handle appropriately. Of course, this is not fool-proof,
% but it should work most of the time.

global SYSSIZE

% Determine the variable names:
str = func2str(inFun);
parLoc = strfind(str, ')');
varList = str(3:parLoc(1)-1);
commaLoc = strfind(varList, ',');
varNames = cell(1, numel(commaLoc)+1);
for k = 1:numel(commaLoc)
    varNames{k} = varList(1:commaLoc(k)-1);
    varList(1:commaLoc(k)) = [];
    commaLoc = commaLoc - commaLoc(k); 
end
varNames{numel(commaLoc)+1} = varList;

% Look for 'diff', 'Diff', or 'D':
diffLoc = cellfun(@(v) any(strcmp(v, {'diff', 'Diff', 'D'})), varNames);
idx = find(diffLoc);
if ( isempty(idx) )
    % None present - new syntax!
    outFun = inFun;
    varNamesOut = varNames;
    return
end

warning('CHEBFUN:CHEBFUN:pde15s:oldSyntax', ...
    ['The syntax\n ' str '\nis deprecated and may not be supported in ', ...
     'future releases.\n Please see the PDE15S documentation for details.'] )

% Switch the order for V5 syntax:
varNamesNewOrder = varNames([(1+SYSSIZE):(idx-1), 1:SYSSIZE]);

% Get the new list of variables:
varList1 = varNamesNewOrder{1};
varList2 = varNames{1};
for k = 2:idx-1
    varList1 = [varList1, ',', varNamesNewOrder{k}]; %#ok<AGROW>
    varList2 = [varList2, ',', varNames{k}]; %#ok<AGROW>
end


% Compile the new function string:
funList = {'@diff', '@sum', '@cumsum', '@fred'};
newStr = ['@(', varList1, ')' 'inFun(' varList2];
for k = 0:(nargin(inFun) - idx)
    newStr = [newStr, ',', funList{k+1}]; %#ok<AGROW>
end
newStr = [newStr, ')'];

% Make the new function handle:
outFun = eval(newStr);

% Return the obtained variable names
varNamesOut = varNamesNewOrder;
end
