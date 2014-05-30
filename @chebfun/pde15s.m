function varargout = pde15s(pdeFun, tt, u0, bc, varargin)
%PDE15S   Solve PDEs using the CHEBFUN system.
%   UU = PDE15s(PDEFUN, TT, U0, BC) where PDEFUN is a handle to a function with
%   arguments u, t, x, and D, TT is a vector, U0 is a CHEBFUN or a CHEBMATRIX,
%   and BC is a chebop boundary condition structure will solve the PDE dUdt =
%   PDEFUN(UU, t, x) with the initial condition U0 and boundary conditions BC
%   over the time interval TT.
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
%   For equations of one variable, UU is output as an array-valued CHEBFUN,
%   where UU(:, k) is the solution at TT(k). For systems, the solution UU is
%   returned as a CHEBMATRIX with the different variables along the rows, and
%   time slices along the columns.
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
%      f = @(u, v, w)  [ .1*diff(u, 2) - 100*u.*v , ...
%                        .2*diff(v, 2) - 100*u.*v , ...
%                        .001*diff(w, 2) + 2*100*u.*v ];
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

% Copyright 2014 by The University of Oxford and The Chebfun Developers. See
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  EVENT SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ctr = 1;
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
        
        for kk = 1:numel(t)
            ctr = ctr + 1;
            % Reshape solution:
            Uk = reshape(U(:,kk), n, SYSSIZE);
            Uk = chebfun(Uk, DOMAIN);
            uOut{ctr} = Uk;

            % Plot solution:
            plotFun(Uk, t(kk));
        end

    end

    function status = adaptiveEvent(t, U, flag)
        % This event is called at the end of each chunk in non-adaptive mode.
        status = false;
        if ( ~isempty(flag) )
            return
        end
        
        for kk = 1:numel(t)
            % Reshape solution:
            Uk = reshape(U(:,kk), currentLength, SYSSIZE);

            Uk2 = Uk*(1:SYSSIZE)'/SYSSIZE;
             f = chebtech2(Uk2, pref2);
%             [ishappy, epslevel, cutoff] = classicCheck(f, Uk2, pref2);

%             C = chebtech2.coeffs2vals(Uk2);
%             C = max(abs(C), [], 2);
%             vscale = norm(U(:), inf);
             [ishappy, epslevel, cutoff] = plateauCheck(f, Uk2, pref2);

            if ( ishappy )  

                % Store these values:
                uCurrent = chebfun(Uk, DOMAIN);
                
                uCoeff = chebpoly(uCurrent);
                first = max(1,size(uCoeff,2)-cutoff);
                uCurrent = chebfun(uCoeff(:,first:end).',DOMAIN,'coeffs');
                uCurrent = simplify(uCurrent,epslevel);
                
                ctr = ctr + 1;
                uOut{ctr} = uCurrent;
                currentT = t;

                % Plot solution:
                plotFun(uCurrent, t(kk));

                if ( cutoff < 0.4*n )
                    currentLength = round(1.25*cutoff)';
                    currentLength = currentLength + 1 - rem(currentLength,2)
                    status = true;
                end

            else 

                % Increase length and bail out:
                currentLength = 2*currentLength-1
                status = true;
            end
        end

    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOTTING SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    function varargout = plotFun(U, t)
        if ( ~doPlot )
            return
        end
        if ( ~guiFlag )
            cla, shg
        end
        set(gcf, 'doublebuf', 'on');

        % Plot:
        h = plot(U, plotOpts{:});

        % Hold?
        ish = ishold();
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
        
        if ( nargin > 1 )
            title(sprintf('%s = %.3f,  len = %i', tlabel, t, length(U)))
        end
        drawnow
        
        if ( nargout > 0 )
            varargout{1} = h;
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ODE tolerances: (AbsTol and RelTol must be <= Tol/10)
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
            A = @(n) [1 zeros(1, n-2) -1]*colloc2.diffmat(n, k)*c;
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
            A = @(n) [1 zeros(1, n-1)]*colloc2.diffmat(n)*diff(DOMAIN)/2;
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
        leftNonlinBCVals = zeros(1,max(sizeOp));  % nominally zero
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
        middleNonlinBCVals = zeros(1,max(sizeOp));  % nominally zero
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
            A = @(n) [zeros(1, n-1) 1]*colloc2.diffmat(n)*diff(DOMAIN)/2;
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
        rightNonlinBCVals = zeros(1,max(sizeOp));  % nominally zero
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
% TODO: These are no longer used?
pref = chebfunpref;
pref.techPrefs.eps = tol;
pref.refinementFunction = 'resampling';
pref.enableBreakpointDetection = 0;
pref.techPrefs.sampleTest = 0;
pref.enableSingularityDetection = 0;
pref2 = chebtech.techPref(pref);
pref2.eps = tol;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% TIME CHUNKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot initial condition:
plotFun(u0, tt(1));

if ( ~isnan(optN) )
    % Non-adaptive in space:
    oneStep(chebpts(optN, DOMAIN), tt); % Do all chunks at once!
    
else
    
    currentLength = max(length(u0), 9);
    currentT = tt(1);
    while ( currentT < tt(end) )
        tt(tt < currentT) = [];
        oneStep(chebpts(currentLength, DOMAIN), tt);
    end
end

if ( doPlot && ~ish )
    hold off
end

% If we only had one dependent variable, return an array valued CHEBFUN instead
% of a QUASIMATRIX.
if ( SYSSIZE == 1 )
    uOut = horzcat(uOut{:});
else
    % TODO: Determine what output we want for systems of equations. CHEBMATRIX?
    blocks = cell(SYSSIZE, numel(uOut));
    for k = 1:SYSSIZE
        blocks(k,:) = cellfun(@(u) extractColumns(u, k), uOut, 'UniformOutput', false);
    end
    uOut = chebmatrix(blocks); % CHEBMATRIX
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

    function U = oneStep(x, tspan)
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
            
            % Projection / mass matrix.
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
        
        % We have to ensure the starting condition satisfies the boundary
        % conditions, or else we will get a singularity in time. We do this by
        % tweaking the BC values to match the reality. Changing the IC itself is
        % much trickier.
        linearBCVals = B*U0(:);
        Utmp = chebdouble(U0, DOMAIN);
        if ( ~isempty(leftNonlinBCLocs) )
            tmp = leftNonlinBCFuns(tspan(1), x, Utmp);
            tmp = double(tmp);
            if ( size(tmp, 1) ~= n )
                tmp = reshape(tmp, n, numel(tmp)/n);
            end
            leftNonlinBCVals = tmp(1, :);
        end
        if ( ~isempty(middleNonlinBCLocs) )
            % TODO: This won't work if there are also left nonlin BCs
            middleNonlinBCVals = double(middleNonlinBCFuns(tspan(1), x, Utmp));
        end
        if ( ~isempty(rightNonlinBCLocs) )
            tmp = rightNonlinBCFuns(tspan(1), x, Utmp);
            tmp = double(tmp);
            if ( size(tmp, 1) ~= n )
                tmp = reshape(tmp, n, numel(tmp)/n);
            end
            rightNonlinBCVals = fliplr(tmp(end, :));
        end       
        
        % ODE options: (mass matrix)
        opt2 = odeset(opt, 'Mass', M, 'MassSingular', 'yes', ...
            'InitialSlope', odeFun(tspan(1), U0), 'MStateDependence', 'none');
        
        % Solve ODE over time chunk with ode15s:
        [ignored, U] = ode15s(@odeFun, tspan, U0, opt2);

        if ( length(tspan) > 2 )
            return
        end
        
        % Reshape solution:
        U = reshape(U(end, :).', n, SYSSIZE);
                
        % The solution we'll take out and store:
        unew = U;
        
        % Collapse systems to single chebfun for constructor. Choose a
        % combination unlikely to lead to cancellation. 
        U = U*(1:SYSSIZE)'/SYSSIZE;
        
        % Zero out irrelevant coefficients, to prevent resolving noise. 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%  ODEFUN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function F = odeFun(t, U)
            % This is what ODE15S() calls.
            %fprintf('  t = %.5f\n',t)
            
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
                    q(l, 1) = feval(BCRHS{l}, t) - linearBCVals(l);
                else
                    q(l, 1) = BCRHS{l} - linearBCVals(l);
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
                F(rows(indx)) = tmp(1, :) - leftNonlinBCVals;
            end
            if ( ~isempty(middleNonlinBCLocs) )
                % TODO: This won't work if there are also left nonlin BCs
                indx = 1:length(middleNonlinBCLocs);
                F(rows(indx)) = double(middleNonlinBCFuns(t, x, Utmp)) - ...
                    middleNonlinBCVals;
            end            
            if ( ~isempty(rightNonlinBCLocs) )
                indx = numel(BCRHS) + 1 - rightNonlinBCLocs;
                tmp = rightNonlinBCFuns(t, x, Utmp);
                tmp = double(tmp);
                if ( size(tmp, 1) ~= n )
                    tmp = reshape(tmp, n, numel(tmp)/n);
                end
                F(rows(indx)) = fliplr(tmp(end, :)) - rightNonlinBCVals;
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

%%

% function [ishappy, epslevel, cutoff] = plateauCheck(coeff, vscale, pref)
% %PLATEAUCHECK   Seek a plateau in Chebyshev coefficients.
% %  Inputs:
% %    coeff:  vector of Chebyshev polynomial coefficients (high order to low)
% %    vscale: indication of the scale to resolve relative to (default=Inf,
% %            no effect)
% %    pref:   cheboppref
% %
% %  Outputs:
% %    ishappy:  true if convergence was achieved
% %    epslevel: the apparent epslevel of the truncation
% %    cutoff:   where to truncate the coefficients
% %
% % This check is needed because of condition numbers in differential equations.
% % We can't be sure that a solution will ever be resolved to full precision, so
% % we have to be willing to stop if the convergence appears to have trailed off.
% 
% % TODO: Unify and locate with the chebtech happiness checks.
% 
% % NaNs are not allowed.
% if ( any(isnan(coeff)) )
%     error('CHEBFUN:FUN:plateauCheck:NaNeval', ...
%         'Function returned NaN when evaluated.')
% end
% 
% % We omit the last 12% because aliasing can pollute them significantly.
% n = length(coeff);
% n88 = ceil( 0.88*n );
% % Preferred tolerance
% epslevel = pref.eps;  
% % Magnitude and rescale.
% if ( vscale > 0 )
%     absCoeff = abs( coeff(n:-1:1) ) / vscale;
% end
% 
% % %%%%%%%%%%%%%%%%%%%%%%%% Serious checking starts here. %%%%%%%%%%%%%%%%%%%%%%%
% % There are two ways to pass the test. Either the coefficients have
% % achieved the goal epslevel, or the convergence appears to have levelled
% % off for good (plateau).
% 
% % Guilty until proven innocent.
% ishappy = false;
% 
% %% 1. Strict test.
% 
% % Find the last place where the coeffs exceed the allowable level.
% % Then go out a bit further to be safe.
% cutoff = 4 + find( absCoeff >= epslevel, 1, 'last' );
% 
% if ( cutoff < 0.95*n88 )
%     % Achieved the strict test.
%     ishappy = true;
%     
% elseif ( n88 < 17 )
%     % If there aren't enough coefficients, give up checking.
%     epslevel = absCoeff(n88);
%     cutoff = n88;
%     
% %% 2. Plateau test.
% else
%     
%     % Demand at least this much accuracy.
%     thresh = max(log(epslevel*100), log(1e-7));
%     
%     % Convergence is usually not far from linear in the log scale.
%     logAbs = log(absCoeff);
%     
%     % Even though some methods can compute really small coefficients relative to
%     % the norm, they ultimately contribute nothing. Also the occasional "zero"
%     % coefficient causes troublesome infinities.
%     logAbs = max( logAbs, log(eps/1000) );
%     
%     % Look for a sustained leveling off in the decrease.
%     
%     % TODO: Use the van Herk filter to do this more efficiently.
%     
%     % Symmetries can cause one or more consecutive coefficients to be zero, and
%     % we only care about the nonzero ones. Use a windowed max to remove the
%     % small values.
%     winSize = 6;
%     winMax = logAbs;
%     for k = 1:winSize
%         winMax = max( winMax(1:end-1), logAbs(k+1:end) );
%     end
%     n88 = length(winMax);
%     
%     %%% Alternative windowed max: This avoids the for loop but might hog memory.
%     %%index = bsxfun(@plus, (1:n)', 0:winsize-1);
%     %%logabs = max(logabs(index),[],2);
%     
%     % Start with a low pass smoothing filter that introduces a lag.
%     lag = 6;
%     LPA = [1, zeros(1,lag-1), -2, zeros(1, lag-1), 1] / (lag^2);
%     LPB = [1, -2, 1];
%     smoothLAC = filter(LPA, LPB, winMax);  % smoothed logabs coeffs
%     
%     % If too little accuracy has been achieved, do nothing.
%     tOK = find(smoothLAC < thresh, 1) - lag;
%     if ( isempty(tOK) || (n88 - tOK < 16) )
%         return
%     end
%     
%     % Smooth the first difference of the smoothed coefficent sequence.
%     smoothDiff = filter( LPA, LPB, diff(smoothLAC) );
%     
%     % Where is the decrease most rapid?
%     SDmin = min(smoothDiff);
%     
%     % Don't look at anything until all substantial decrease has ended.
%     tstart = find( smoothDiff < 0.25*SDmin, 1, 'last' );
%     
%     % Find where the decrease has permanently slowed to 10% of the fastest.
%     isSlow = smoothDiff(tstart:end) > 0.01*SDmin;
%     lastFast = find(~isSlow, 1, 'last');
%     if ( isempty(lastFast) )
%         lastFast = 0;
%     end
%     slow = tstart + lastFast - 1 + find( isSlow(lastFast+1:end) );
%     slow = slow - floor(lag/2);  % compensate for the filter lag
%     
%     % Find the first run of 5 consecutive slow hits.
%     first = find( slow(5:end) - slow(1:end-4) == 4, 1 );  % may be empty
%     cutoff = slow(first);  % may be empty, will give false next
%     
%     % If the cut location is within the coefficient sequence, we're done.
%     if ( cutoff < n88 )
%         ishappy = true;
%     end
%     
% end
% 
% if ( ishappy )
%     % Use the information from the cut to deduce an eps level.
%     winEnd = min( n88, cutoff + 4 );
%     epslevel = max( absCoeff(cutoff:winEnd) );
% end
% 
% end
