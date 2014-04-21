function varargout = pde15s(pdeFun, tt, u0, bc, varargin)
%PDE15S   Solve PDEs using the CHEBFUN system.
%   UU = PDE15s(PDEFUN, TT, U0, BC) where PDEFUN is a handle to a function with
%   arguments u, t, x, and D, TT is a vector, U0 is a chebfun, and BC is a
%   chebop boundary condition structure will solve the PDE dUdt = PDEFUN(UU, t,
%   x) with the initial condition U0 and boundary conditions BC over the time
%   interval TT.
%
%   PDEFUN should take the form @(T, X, U1, U2, ..., UN), where U1, ..., UN are
%   the unknown dependent variables to be solved for, T is time, and X is space.
%
%   For backwards compatability, the syntax @(U1, U2, ..., UN, T, X, D, S, C)
%   for PDEFUN, where U1, ..., UN are the unknown dependent variables to be
%   solved for, T is time, X is space, D is the differential operator, S is the
%   definite integral operator (i.e., 'sum'), and C the indefinite integral
%   operator (i.e., 'cumsum') is also supported.
%
%   For equations of one variable, UU is output as an array valued chebfun,
%   where UU(:, k) is the solution at TT(k). For systems, the solution is
%   returned as a cell array of array valued chebfuns.
%
% Example 1: Nonuniform advection
%     x = chebfun('x', [-1 1]);
%     u = exp(3*sin(pi*x));
%     f = @(t, x, u) -(1 + 0.6*sin(pi*x)).*diff(u) + 5e-5*diff(u, 2);
%     opts = pdeset('Ylim', [0 20], 'PlotStyle', {'LineWidth', 2});
%     uu = pde15s(f, 0:.05:3, u, 'periodic', opts);
%     surf(uu, 0:.05:3)
%
% Example 2: Kuramoto-Sivashinsky
%     x = chebfun('x');
%     u = 1 + 0.5*exp(-40*x.^2);
%     bc.left = @(u) [u - 1, diff(u)];
%     bc.right = @(u) [u - 1, diff(u)];
%     f = @(u) u.*diff(u) - diff(u, 2) - 0.006*diff(u, 4);
%     opts = pdeset('Ylim', [-30 30], 'PlotStyle', {'LineWidth', 2});
%     uu = pde15s(f, 0:.01:.5, u, bc, opts);
%     surf(uu, 0:.01:.5)
%
% Example 3: Chemical reaction (system)
%      x = chebfun('x');
%      u = [ 1 - erf(10*(x+0.7)) , 1 + erf(10*(x-0.7)) , 0 ];
%      f = @(u, v, w, diff)  [ .1*diff(u, 2) - 100*u.*v , ...
%                           .2*diff(v, 2) - 100*u.*v , ...
%                           .001*diff(w, 2) + 2*100*u.*v ];
%      opts = pdeset('Ylim', [0 2], 'PlotStyle', {'LineWidth', 2});
%      uu = pde15s(f, 0:.1:3, u, 'neumann', opts);
%      mesh(uu{3})
%
% See chebfun/test/test_pde15s.m for more examples.
%
%   UU = PDE15s(PDEFUN, TT, U0, BC, OPTS) will use nondefault options as defined
%   by the structure returned from OPTS = PDESET.
%
%   UU = PDE15s(PDEFUN, TT, U0, BC, OPTS, N) will not adapt the grid size in
%   space. Alternatively OPTS.N can be set to the desired size.
%
%   [TT, UU] = PDE15s(...) returns also the time chunks TT.
%
%   There is some support for nonlinear and time-dependent boundary conditions,
%   such as
%       x = chebfun('x', [-1 1]);
%       u = exp(-3*x.^2);
%       f = @(t, x, u) .1*diff(u, 2);
%       bc.left = @(t, x, u) u - t;
%       bc.right = 0;
%       opts = pdeset('Ylim', [0 2], 'PlotStyle', {'LineWidth', 2});
%       uu = pde15s(f, 0:.1:2, u, bc, opts);
%       waterfall(u);
%   with the input format being the same as PDEFUN described above.
%
% See also PDESET, ODE15S.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% TODO: Syncronise with CHEBOP syntax. (In particular, .lbc, .rbc, and .bc).

global DIFFORDER SYSSIZE
DIFFORDER = 0;
SYSSIZE = 0;

% Default options:
tol = 1e-6;             % 'eps' in Chebfun terminology
doPlot = 1;             % Plot after every time chunk?
doHold = 0;             % Hold plot?
plotOpts = {'-'};       % Plotting style

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

% PDE solver options:
if ( ~isempty(opt.Eps) )
    tol = opt.Eps;
end
if ( ~isempty(opt.Plot) )
    doPlot = strcmpi(opt.Plot, 'on');
end
if ( ~isempty(opt.HoldPlot) )
    doHold = strcmpi(opt.HoldPlot, 'on');
end
if ( ~isempty(opt.PlotStyle) )
    plotOpts = opt.PlotStyle;
end

% Experimental feature for coupled ode/pde systems: (An entry equal to 1 denotes
% that the corresponding variable appears with a time derivative. 0 otherwise.)
if ( isfield(opt, 'PDEflag') )
    pdeFlag = opt.PDEflag;
else
    pdeFlag = true;
end

% Determine which figure to plot to (for CHEBGUI) and set default display values
% for variables.
% TODO: Ensure this still works when CHEBGUI gets merged into development.
YLim = opt.YLim;
gridOn = 0;
guiFlag = false;
if ( isfield(opt, 'handles') )
    if ( opt.handles.gui )
        guiFlag = true;
        axesSol = opt.handles.fig_sol;
        axesNorm = opt.handles.fig_norm;
        axes(axesSol);
        gridOn = opt.handles.guifile.options.grid;
        solveButton = opt.handles.button_solve;
        clearButton = opt.handles.button_clear;
    end
    varNames = opt.handles.varnames;
    xLabel = opt.handles.indVarName{1};
    tlabel = opt.handles.indVarName{2};
else
    varNames = 'u';
    xLabel = 'x';
    tlabel = 't';
end

% ODE tolerances: (AbsTol and RelTol must be <= Tol/10)
aTol = odeget(opt, 'AbsTol', tol/10);
rTol = odeget(opt, 'RelTol', tol/10);
if ( isnan(optN) )
    aTol = min(aTol, tol/10);
    rTol = min(rTol, tol/10);
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
            error('CHEBFUN:pde15s:piecewise', ...
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

% Get the domain and the independent variable 'x':
DOMAIN = domain(u0, 'ends');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%  PARSE INPUTS TO PDEFUN  %%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the size of the system, i.e. number of dependent variables.
SYSSIZE = min(size(u0));
pdeFun = parseFun(pdeFun);
if ( isfield(opt, 'difforder') )
    DIFFORDER = opt.difforder;
else
    getDIFFORDER(pdeFun);
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

if ( ischar(bc) && strcmpi(bc, 'periodic') )
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% PERIODIC BCS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r = cell(sum(DIFFORDER), 1);
    count = 1;
    for j = 1:SYSSIZE
        for k = 1:DIFFORDER(j)
            c = (diff(DOMAIN)/2)^k;
            A = @(n) [1 zeros(1, n-2) -1]*diffmat(n, k)*c;
            r{count} = @(n) [zeros(1, (j-1)*n) A(n) zeros(1,(SYSSIZE-j)*n)];
            count = count + 1;
        end
    end
    bc = struct( 'left', [], 'right', []);
    bc.left.op = r(1:2:end);
    bc.right.op = r(2:2:end);
    BCRHS = num2cell(zeros(1, numel(r)));
    
else
    
    if ( isfield(bc, 'left') && ~isfield(bc, 'right') )
        bc.right = [];
    elseif ( isfield(bc, 'right') && ~isfield(bc, 'left') )
        bc.left = [];
    elseif ( ~isfield(bc, 'left') && ~isfield(bc, 'right') )
        bc = struct('left', [], 'right', [], 'middle', bc);
%         bc.left = bc;
%         bc.right = bc.left;
    end
    
    % Deal with struct and numeric input:
    bc.left = dealWithStructInput(bc.left);
    bc.right = dealWithStructInput(bc.right);
    
    if ( isempty(bc.left) )
        bc.left.op = [];
    elseif ( ischar(bc.left) || (iscell(bc.left) && ischar(bc.left{1})) )
        %% %%%%%%%%%%%%%%%%%%%%% DIRICHLET AND NEUMANN BCS (LEFT) %%%%%%%%%%%%%%%%%%
        if ( iscell(bc.left) )
            v = bc.left{2};
            bc.left = bc.left{1};
        else
            v = 0;
        end
        if ( ~isnumeric(v) )
            error('For BCs of the form {char, val} val must be numeric.')
        end
        if ( strcmpi(bc.left, 'dirichlet') )
            A = @(n) [1 zeros(1, n - 1)];
        elseif ( strcmpi(bc.left, 'neumann') )
            % TODO: Make left diff operator explicitly.
            A = @(n) [1 zeros(1, n-1)]*diffmat(n)*diff(DOMAIN)/2;
        else
            error('Unknown BC syntax');
        end
        bc.left.op = cell(SYSSIZE, 1);
        for k = 1:SYSSIZE
            bc.left.op{k} = @(n) [zeros(1, ( k -1)*n) A(n) ...
                zeros(1 , (SYSSIZE - k)*n)];
        end
        BCRHS = num2cell(repmat(v, SYSSIZE, 1));
    elseif ( numel(bc.left) == 1 && isa(bc.left, 'function_handle') )
        %% %%%%%%%%%%%%%%%%%%%%%%%%%  GENERAL BCS (LEFT)  %%%%%%%%%%%%%%%%%%%%%%%%%%
        op = parseFun(bc.left);
        tmp = chebdouble(ones(1, SYSSIZE));
        sizeOp = size(op(0, mean(DOMAIN), tmp));
        leftNonlinBCLocs = 1:max(sizeOp);
        bc.left = [];
        bc.left.op = {@(n) zeros(max(sizeOp), SYSSIZE*n)}; % Dummy entries.
        BCRHS = num2cell(zeros(1, max(sizeOp)));
        leftNonlinBCFuns = op;
    else
        error('Unknown BC syntax');
    end
    
    if ( isfield(bc, 'middle') && isa(bc.middle, 'function_handle') )
        %% %%%%%%%%%%%%%%%%%%%%%  GENERAL BCS (MIDDLE)  %%%%%%%%%%%%%%%%%%%%%%%%
        op = parseFun(bc.middle);
        tmp = chebdouble(ones(1, SYSSIZE));
        sizeOp = size(op(0, mean(DOMAIN), tmp));
        middleNonlinBCLocs = 1:max(sizeOp);
        bc.middle = [];
        bc.middle.op = {@(n) zeros(max(sizeOp), SYSSIZE*n)}; % Dummy entries.
        BCRHS = num2cell(zeros(1, max(sizeOp)));
        middleNonlinBCFuns = op;
    end
    
    if ( isempty(bc.right) )
        bc.right.op = [];
    elseif ( ischar(bc.right) || (iscell(bc.right) && ischar(bc.right{1})) )
        %% %%%%%%%%%%%%%%%%%%%%% DIRICHLET AND NEUMANN BCS (RIGHT) %%%%%%%%%%%%%%%%%
        if ( iscell(bc.right) )
            v = bc.right{2};
            bc.right = bc.right{1};
        else
            v = 0;
        end
        if ( ~isnumeric(v) )
            error('For BCs of the form {char, val} val must be numeric.')
        end
        
            
        if ( strcmpi(bc.right, 'dirichlet') )
            A = @(n) [zeros(1, n-1), 1];
        elseif ( strcmpi(bc.right, 'neumann') )
            % TODO: Make right diff operator explicitly.
            A = @(n) [zeros(1, n-1) 1]*diffmat(n)*diff(DOMAIN)/2;
        else
            error('Unknown BC syntax');
        end
        bc.right.op = cell(SYSSIZE, 1);
        for k = 1:SYSSIZE
            bc.right.op{k} = @(n) [zeros(1,(k-1)*n) A(n) zeros(1,(SYSSIZE-k)*n)];
        end
        BCRHS = [BCRHS num2cell(repmat(v, SYSSIZE, 1))];
        
    elseif ( numel(bc.right) == 1 && isa(bc.right, 'function_handle') )
        %% %%%%%%%%%%%%%%%%%%%%%%%%  GENERAL BCS (RIGHT)  %%%%%%%%%%%%%%%%%%%%%%%%%%
        op = parseFun(bc.right);
        tmp = chebdouble(ones(1, SYSSIZE));
        sizeOp = size(op(0, mean(DOMAIN), tmp));
        rightNonlinBCLocs = 1:max(sizeOp);
        bc.right = [];
        bc.right.op = {@(n) zeros(max(sizeOp), SYSSIZE*n)};
        BCRHS = [BCRHS num2cell(zeros(1, max(sizeOp)))];
        rightNonlinBCFuns = op;
    else
        error('Unknown BC syntax');
    end
    
end

if ( ~isfield(bc, 'middle') )
    bc.middle.op = [];
end

%% %%%%%%%%%%%%%%%%%%%%%%%% Support for coupled BVP-PDEs! %%%%%%%%%%%%%%%%%%%%%%
% TODO: Test this!
if ( ~all(pdeFlag) )
    userMassSet = true;
    userMass = @(n) [];
    for k = 1:numel(pdeFlag)
        if ( pdeFlag(k) )
            A = @(n) eye(n);
        else
            A = @(n) zeros(n);
        end
        userMass = @(n) [userMass(n) ; ...
            zeros(n,(k-1)*n), A(n), zeros(n,(SYSSIZE-k)*n)];
    end
else
    userMassSet = false;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOTTING SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( doPlot )
    
    if ( ~guiFlag )
        cla, shg
    end
    set(gcf, 'doublebuf', 'on');
    
    % Plot:
    plot(u0, plotOpts{:});
    
    % Hold?
    ish = ishold(); % Store to be able to return to previous state once finished
    if ( doHold )
        hold on
    end
    
    % Fix Y limits?
    if ( ~isempty(YLim) )
        ylim(YLim);
    end
    
    % Axis labels and legends:
    xlabel(xLabel);
    if ( numel(varNames) > 1 )
        legend(varNames);
    else
        ylabel(varNames);
    end
    
    % Grid?
    if ( gridOn )
        grid on
    end
    drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The vertical scale of the intial condition:
vscl = get(u0, 'vscale');

% Initial condition:
uCurrent = u0;
% Storage:
uOut = cell(1, numel(tt));
uOut{1} = uCurrent;

% Initialise variables for ONESTEP():
B = []; q = []; rows = []; M = []; n = [];

% Set the preferences:
pref = chebpref;
pref.techPrefs.eps = tol;
pref.refinementFunction = 'resampling';
pref.enableBreakpointDetection = 0;
pref.techPrefs.sampleTest = 0;
pref.enableSingularityDetection = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% TIME CHUNKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin time chunks
for nt = 1:length(tt)-1
    
    % Solve one chunk:
    if ( isnan(optN) )
        % Size of current length:
        currentVscale = vscale(uCurrent);
        tmp = max(currentVscale, vscl) - currentVscale;
        currentLength = length(simplify(uCurrent + tmp, tol));
        pref.techPrefs.minPoints = max(2*currentLength, 9);
        vscl = max([vscl, currentVscale]);
        chebfun( @(x) vscl + oneStep(x), DOMAIN, pref);
    else
        % Non-adaptive in space:
        currentLength = optN;
        oneStep(chebpts(optN, DOMAIN));
    end
    
    % Get chebfun of solution from this time chunk:
    uCurrent = chebfun(unew, DOMAIN);
    
    % Store in uOut:
    uOut{nt + 1} = uCurrent;
    
    % Plotting:
    if ( doPlot )
        plot(uCurrent, plotOpts{:});
        if ( ~isempty(YLim) )
            ylim(YLim);
        end
        if ( ~doHold )
            hold off
        end
        % Axis labels
        xlabel(xLabel);
        if ( numel(varNames) > 1 )
            legend(varNames);
        else
            ylabel(varNames);
        end
        % Determines whether grid is on
        if ( gridOn )
            grid on
        end
        title(sprintf('%s = %.3f,  len = %i', tlabel, tt(nt+1), currentLength)), drawnow
        %     elseif ( guiFlag )
        %         drawnow
    end
    
    % TODO: Restore once CHEBGUI is operating.
    %     if ( guiFlag )
    %         % Interupt comutation if stop or pause  button is pressed in the GUI.
    %         if ( strcmp(get(solveButton, 'String'), 'Solve') )
    %             tt = tt(1:nt+1);
    %             if SYSSIZE == 1,
    %                 uOut = uOut(1:nt+1);
    %             else
    %                 for k = 1:SYSSIZE
    %                     uOut{k} = uOut{k}(1:nt+1);
    %                 end
    %             end
    %             break
    %         elseif ( strcmp(get(clearButton, 'String'), 'Continue') )
    %             defaultlinewidth = 2;
    %             axes(axesNorm)
    %             if ( ~iscell(uOut) )
    %                 waterfall(uOut(1:nt+1), tt(1:nt+1), 'simple', 'linewidth', defaultlinewidth)
    %                 xlabel(xLabel), ylabel(tlabel), zlabel(varnames)
    %             else
    %                 cols = get(0, 'DefaultAxesColorOrder');
    %                 for k = 1:numel(uOut)
    %                     plot(0, NaN, 'linewidth', defaultlinewidth, 'color', cols(k, :)), hold on
    %                 end
    %                 legend(varnames);
    %                 for k = 1:numel(uOut)
    %                     waterfall(uOut{k}, tt(1:nt+1), 'simple', 'linewidth', ...
    %                           defaultlinewidth, 'edgecolor', cols(k, :)), hold on
    %                     xlabel(xLabel), ylabel(tlabel)
    %                 end
    %                 view([322.5 30]), box off, grid on, hold off
    %             end
    %             axes(axesSol)
    %             waitfor(clearButton, 'String');
    %         end
    %     end
end

if ( doPlot && ~ish )
    hold off
end

% If we only had one dependent variable, return an array valued CHEBFUN instead
% of a QUASIMATRIX.
% TODO: Determine what we want to output for systems of equations. CHEBMATRIX?
if ( SYSSIZE == 1 )
    uOut = horzcat(uOut{:});
end

switch nargout
    case 0
    case 1
        varargout{1} = uOut;
    case 2
        varargout{1} = tt;
        varargout{2} = uOut;
    otherwise
        varargout{1} = tt;
        varargout{2} = uOut;
        varargout{3:nargout} = [];
        warning('CHEBFUN:pde15s:output', ...
            'PDE15S has a maximum of two outputs.');
end

clear global DIFFORDER
clear global SYSSIZE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ONESTEP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function U = oneStep(x)
        % Constructs the result of one time chunk at fixed discretization.
        
        if ( length(x) == 2 )
            U = [0 ; 0];
            return
        end
        
        % Evaluate the chebfun at discrete points:
        U0 = feval(uCurrent, x);
        
        % This depends only on the size of n. If this is the same, reuse!
        if ( isempty(n) || n ~= length(x) )
            
            % The new discretisation length
            n = length(x);
            
            % Linear constraints:
            bcop = [bc.left.op ; bc.middle.op ; bc.right.op];
            B = cell2mat(cellfun(@(f) feval(f, n), bcop, 'UniformOutput', false));
            
            % Project / mass matrix.
            M = {};
            for kk = 1:SYSSIZE
                xk = chebpts(n-DIFFORDER(kk), DOMAIN, 1);
                M{kk} = barymat(xk, x);
            end
            M = [ 0*B ; blkdiag(M{:})];
            rows = 1:size(B, 1);
            
            % Multiply by user-defined mass matrix
            if ( userMassSet )
                M = feval(userMass, n)*M;
            end
            
        end
        
        % ODE options: (mass matrix)
        opt2 = odeset(opt, 'Mass', M, 'MassSingular', 'yes', ...
            'InitialSlope', odeFun(tt(nt), U0), 'MStateDependence', 'none');
        
        % Solve ODE over time chunk with ode15s:
        [ignored, U] = ode15s(@odeFun, tt(nt:nt+1), U0, opt2);
        
        % Reshape solution:
        U = reshape(U(end, :).', n, SYSSIZE);
        
        % The solution we'll take out and store:
        unew = U;
        
        % Collapse systems to single chebfun for constructor (is addition right?)
        U = sum(U, 2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%  ODEFUN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function F = odeFun(t, U)
            % This is what ODE15S() calls.
            
            % Reshape to n by SYSSIZE:
            U = reshape(U, n, SYSSIZE);
            
            % Evaluate the PDEFUN:
            Utmp = chebdouble(U, DOMAIN);
            F = pdeFun(t, x, Utmp);
            F = double(F);
            F = M*F(:);
            
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
                tmp = leftNonlinBCFuns(t, x, Utmp);
                tmp = double(tmp);
                if ( size(tmp, 1) ~= n )
                    tmp = reshape(tmp, n, numel(tmp)/n);
                end
                F(rows(indx)) = tmp(1, :);
            end
            if ( ~isempty(middleNonlinBCLocs) )
                % TODO: This won't work if there are also left nonlin BCs
                indx = 1:length(middleNonlinBCLocs);
                F(rows(indx)) = double(middleNonlinBCFuns(t, x, Utmp));
            end            
            if ( ~isempty(rightNonlinBCLocs) )
                indx = numel(BCRHS) + 1 - rightNonlinBCLocs;
                tmp = rightNonlinBCFuns(t, x, Utmp);
                tmp = double(tmp);
                if ( size(tmp, 1) ~= n )
                    tmp = reshape(tmp, n, numel(tmp)/n);
                end
                F(rows(indx)) = fliplr(tmp(end, :));
            end
            
            % Reshape to back to a single column:
            F = F(:);
        end
    end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%  MISC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outFun = parseFun(inFun)
% Rewrites the input function handle to call the right DIFF, SUM, methods, etc,
% and convert the input @(t, x, u, v, w, ...) to @(t, x, U).

global SYSSIZE

% Ensure backwards compatability by procesing the input function.
inFun = backCompat(inFun);

% V4: We don't accept only time or space as input args (i.e., both or nothing).
% V5: Actually, we now accept @(u, x) diff(u, 2) + sin(x).*u (for chebops).
Nind = nargin(inFun) - SYSSIZE;
if ( Nind == 0 )
    outFun = @(t, x, u) conv2cell(inFun, u);
elseif ( Nind == 1 )
    outFun = @(t, x, u) conv2cell(inFun, u, x);
elseif ( Nind == 2 )
    outFun = @(t, x, u) conv2cell(inFun, u, t, x);
else
    error('CHEBFUN:pde15s:inputs_ind', ['Incorrect number of independant ' ...
        'variables in input function. (Must be 0, 1, or 2).']);
end

end

function newFun = conv2cell(oldFun, u, varargin)
% This function allows the use of different variables in the anonymous
% function, rather than using the clunky quasi-matrix notation.

global SYSSIZE

% Note. This is faster than mat2cell!
tmpCell = cell(1, SYSSIZE);

d = get(u, 'domain');
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
if ( isstruct(in) )
    if ( numel(in) == 1)
        warning('CHEBFUN:pde15s:bcstruct', ...
            'PDE15S no longer supports struct inputs for bc.left and bc.right.')
        if ( isfield(in, 'val') )
            out = {in.op, in.val};
        else
            out = in.op;
        end
    else
        error('CHEBFUN:pde15s:bcstruct', ...
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

function outFun = backCompat(inFun)
% In V4 PDE15S required the user to pass in dummy function handles for the
% differential operators. For example, u'' + u' + sum(u) would have needed
%  pdefun = @(u, t, x, diff, sum) diff(u, 2) + diff(u) + sum(u).
% V5 deals with this in a different way (by using the chebdouble class), but we
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
    return
end

warning('CHEBFUN:pde15s:oldSyntax', ...
    ['The syntax\n ' str '\nis depricated and may not be supported in future releases.\n', ...
    'Please see the PDE15S documentation for details.'] )

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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: DIFFMAT and CUMSUMMAT are in @colloc2/private/. Below should be removed.

function D = diffmat(N,k)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
% DIFFMAT  Chebyshev differentiation matrix
% D = DIFFMAT(N) is the matrix that maps function values at N Chebyshev
% points to values of the derivative of the interpolating polynomial at
% those points.
%
% D = DIFFMAT(N,K) is the same, but for the Kth derivative.
%
% The matrices are computed using the 'hybrid' formula of Schneider &
% Werner [1] and Welfert [2] proposed by Tee [3].

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% References:
%  [1] Schneider, C. and Werner, W., "Some new aspects of rational
%   interpolation", Math. Comp. (47) 285--299, 1986.
%  [2] Welfert, B. D., "Generation of pseudospectral matrices I", SINUM,
%   (34) 1640--1657.
%  [3] Tee, T. W., "An adaptive rational spectral method for differential
%   equations with rapidly varying solutions", Oxford DPhil Thesis, 2006.

if nargin < 2, k = 1; end
if N == 0, D = []; return, end
if N == 1, D = 0; return, end

% construct Chebyshev grid and weights
x = chebtech2.chebpts(N);
w = [.5 ; ones(N-1,1)]; w(2:2:end) = -1; w(N) = .5*w(N);

ii = (1:N+1:N^2)';              % indices of diagonal
Dx = bsxfun(@minus,x,x');       % all pairwise differences
Dx(ii) = Dx(ii) + 1;            % add identity
Dxi = 1./Dx;                    % reciprocal
Dw = bsxfun(@rdivide,w.',w);    % pairwise divisions
Dw(ii) = Dw(ii) - 1;            % subtract identity

% k = 1
D = Dw .* Dxi;
D(ii) = 0; D(ii) = - sum(D,2);              % negative sum trick

if k == 1, return, end

% k = 2
D = 2*D .* (repmat(D(ii),1,N) - Dxi);
D(ii) = 0; D(ii) = - sum(D,2);              % negative sum trick

% higher orders
for n = 3:k
    D = n*Dxi .* (Dw.*repmat(D(ii),1,N) - D);
    D(ii) = 0; D(ii) = - sum(D,2);          % negative sum trick
end

end
