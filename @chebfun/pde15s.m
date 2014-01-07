function varargout = pde15s(pdeFun, tt, u0, bc, varargin)
%PDE15S   Solve PDEs using the CHEBFUN system.
%   UU = PDE15s(PDEFUN, TT, U0, BC) where PDEFUN is a handle to a function with
%   arguments u, t, x, and D, TT is a vector, U0 is a chebfun, and BC is a
%   chebop boundary condition structure will solve the PDE dUdt = PDEFUN(UU, t,
%   x) with the initial condition U0 and boundary conditions BC over the time
%   interval TT.
%
%   PDEFUN should take the form @(U1, U2, ..., UN, T, X, D, S, C), where U1,
%   ..., UN are the unknown dependent variables to be solved for, T is time, X
%   is space, D is the differential operator, S is the definite integral
%   operator (i.e., 'sum'), and C the indefinite integral operator (i.e.,
%   'cumsum').
%
%   For equations of one variable, UU is output as a quasimatrix, where UU(:, k)
%   is the solution at TT(k). For systems, the solution is returned as a cell
%   array of quasimatrices.
%
% Example 1: Nonuniform advection
%     x = chebfun('x', [-1 1]);
%     u = exp(3*sin(pi*x));
%     f = @(u, t, x, diff) -(1+0.6*sin(pi*x)).*diff(u) + 5e-5*diff(u,2);
%     opts = pdeset('Ylim', [0 20], 'PlotStyle', {'LineWidth', 2});
%     uu = pde15s(f, 0:.05:3, u, 'periodic', opts);
%     surf(uu, 0:.05:3)
%
% Example 2: Kuramoto-Sivashinsky
%     x = chebfun('x');
%     u = 1 + 0.5*exp(-40*x.^2);
%     bc.left = @(u, diff) [u - 1, diff(u)];
%     bc.right = @(u, diff) [u - 1, diff(u)];
%     f = @(u, diff) u.*diff(u) - diff(u, 2) - 0.006*diff(u, 4);
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
%       f = @(u, t, x, diff) .1*diff(u, 2);
%       bc.left = @(u, t, x, diff) u - t;
%       bc.right = 0;
%       opts = pdeset('Ylim', [0 2], 'PlotStyle', {'LineWidth', 2});
%       uu = pde15s(f, 0:.1:2, u, bc, opts);
%       waterfall(u);
%   with the input format being the same as PDEFUN described above.
%
% See also PDESET, ODE15S.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

global DIFFORDER GLOBX DMAT DOMAIN SYSSIZE
DIFFORDER = 0;
GLOBX = [];
DMAT = {};
DOMAIN = [];
SYSSIZE = 0;

% Default options:
tol = 1e-6;             % 'eps' in chebfun terminology
doPlot = 1;             % plot after every time chunk?
doHold = 0;             % hold plot?
plotOpts = {'-'};         % Plot Style

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

% Experimental feature for coupled ode/pde systems:
if ( isfield(opt, 'PDEflag') )
    pdeFlag = opt.PDEflag;
else
    pdeFlag = true;
end

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
DOMAIN = domain(u0);
xd = chebfun(@(x) x, DOMAIN);

% These are used often:
Z = linop.zeros(DOMAIN);
z = linop.zero(DOMAIN);
I = linop.eye(DOMAIN);
D = linop.diff(DOMAIN);
Eleft = linop.feval(DOMAIN(1), DOMAIN);
Eright = linop.feval(DOMAIN(end), DOMAIN);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%  PARSE INPUTS TO PDEFUN  %%%%%%%%%%%%%%%%%%%%%%%%%

% Determine the size of the system
SYSSIZE = min(size(u0));  
pdeFun = parseFun(pdeFun);
getDIFFORDER(pdeFun)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%% PARSE BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( ischar(bc) && (strcmpi(bc, 'neumann') || strcmpi(bc, 'dirichlet')) )
    bc = struct( 'left', bc, 'right', bc);
elseif ( iscell(bc) && numel(bc) == 2 )
    bc = struct( 'left', bc{1}, 'right', bc{2});
end

% Initialise some rubbish:
nllbc = []; nlbcs = {}; GLOBX = 1; funFlagL = false; rhs = {};
nlrbc = []; numlbc = 0; funFlagR = false;

if ( ischar(bc) && strcmpi(bc, 'periodic') )
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%% PERIODIC BCS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r = {};
    count = 1;
    for j = 1:SYSSIZE
        for k = 1:DIFFORDER(j)
            Dk = linop.diff(DOMAIN, k);
            A = Eleft*Dk - Eright*Dk;
            r{count} = [repmat(z, 1, j-1) A repmat(z, 1, SYSSIZE-j)];
            count = count + 1;
        end
    end
    bc = struct( 'left', [], 'right', []);
    bc.left.op = vertcat(r{1:2:end});
    bc.right.op = vertcat(r{2:2:end});
    rhs = num2cell(zeros(1, numel(r)));
    
    
else

    if ( isfield(bc, 'left') && ~isfield(bc, 'right') )
        bc.right = [];
    elseif ( isfield(bc, 'right') && ~isfield(bc, 'left') )
        bc.left = [];
    end
    
    bc.left = dealWithStructInput(bc.left);
    bc.right = dealWithStructInput(bc.right);
    

    %% %%%%%%%%%%%%%%%%%%%%% DIRICHLET AND NEUMANN BCS  %%%%%%%%%%%%%%%%%%%%%%%%
    
    if ( ischar(bc.left) || (iscell(bc.left) && ischar(bc.left{1})) )
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
            A = Eleft;
        elseif ( strcmpi(bc.left, 'neumann') )
            A = Eleft*D;
        end
        op = cell(SYSSIZE, 1);
        for k = 1:SYSSIZE
            op{k} = [repmat(z, 1, k-1) A repmat(z, 1, SYSSIZE-k)];
        end
        bc.left.op = vertcat(op{:});
        rhs = num2cell(repmat(v, SYSSIZE, 1));
    end
    if ( ischar(bc.right) || (iscell(bc.right) && ischar(bc.right{1})) )
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
            A = Eright;
        elseif ( strcmpi(bc.right, 'neumann') )
            A = Eright*D;
        end
        op = cell(SYSSIZE, 1);
        for k = 1:SYSSIZE
            op{k} = [repmat(z, 1, k-1) A repmat(z, 1, SYSSIZE-k)];
        end
        bc.right.op = vertcat(op{:});
        rhsTmp = num2cell(repmat(v, SYSSIZE, 1)); % Used in case 0 of gen bs(r).
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%  GENERAL BCS (LEFT)  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ( isfield(bc.left, 'op') && isa(bc.left.op, 'chebmatrix') )
        % 0) Do nothing.
    elseif ( numel(bc.left) == 1 && isa(bc.left, 'function_handle') )
        % 1) Deal with the case where bc is a function_handle vector:

        op = parseFun(bc.left, 'flag');
        sizeOp = size(op(ones(1, SYSSIZE), 0, mean(DOMAIN)));
        nllbc = 1:max(sizeOp);
        bc.left = struct( 'op', []);

        % Dummy entries:
        rowsl = cell(max(sizeOp), 1);
        for k = nllbc
            if ( SYSSIZE == 1 )
                rowsl{k} = Eleft;
            else
                rowsl{k} = [repmat(z, 1, k-1) Eleft repmat(z, 1, SYSSIZE-k)];
            end
        end
        bc.left.op = vertcat(rowsl{:});

        rhs = num2cell(zeros(1, max(sizeOp)));
        nlbcsl = op;
        funFlagL = true;

    elseif ( numel(bc.left) > 0 )
        % 2) Deal with other forms of input

        if ( isa(bc.left, 'linop') )
            rhs = num2cell(zeros(1, size(bc.left,1)));
            bc.left = struct( 'op', bc.left);
        elseif ( isnumeric(bc.left) )
            rhs = num2cell(bc.left);
            bc.left = struct( 'op', Eleft);
            
        else
            error('Unkown BC type');
        end
        
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%  GENERAL BCS (RIGHT)  %%%%%%%%%%%%%%%%%%%%%%%%%%
    numlbc = numel(rhs);
    if ( isfield(bc.right, 'op') && isa(bc.right.op, 'chebmatrix') )
        % 0) Do nothing
        rhs = [rhs rhsTmp];
    elseif ( numel(bc.right) == 1 && isa(bc.right, 'function_handle') )
        % 1) Deal with the case where bc is a function handle vector
        op = parseFun(bc.right, 'flag');
        sizeOp = size(op(ones(1, SYSSIZE), 0, mean(DOMAIN)));
        nlrbc = 1:max(sizeOp);
        bc.right = struct( 'op', []);

        % Dummy entries:
        E = linop.feval(DOMAIN(end), DOMAIN);
        rowsr = cell(max(sizeOp), 1);
        for k = nlrbc
            if ( SYSSIZE == 1 )
                rowsr{k} = E;
            else
                rowsr{k} = [repmat(z, 1, k-1) E repmat(z, 1, SYSSIZE-k)];
            end
        end
        bc.right.op = vertcat(rowsr{:});
        bc.right.val = zeros(max(sizeOp),1);

        rhs = [rhs num2cell(zeros(1, max(sizeOp)))];
        nlbcsr = op;
        funFlagR = true;

    elseif ( numel(bc.right) > 0 )
        % 2) Deal with other forms of input

        if ( isa(bc.right, 'linop') )
            rhs = num2cell(zeros(1, size(bc.right,1)));
            bc.right = struct( 'op', bc.right);
        elseif ( isnumeric(bc.right) )
            rhs = num2cell(bc.right);
            bc.right = struct( 'op', Eright);
        end

    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%% Support for coupled BVP-PDEs! %%%%%%%%%%%%%%%%%%%%%%
if ( ~all(pdeFlag) )
    userMassSet = true;
    userMass = [];
    for k = 1:numel(pdeFlag)
        if ( pdeFlag(k) )
            A = I;
        else
            A = Z;
        end
        userMass = [userMass ; repmat(Z, 1, k-1) A repmat(Z, 1, SYSSIZE-k)];
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
    ish = ishold();
    if ( doHold )
        hold on
    end
    
    % Y limits?
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
% This is needed inside the nested function onestep()
diffOp = repmat(Z, SYSSIZE, SYSSIZE);
for k = 1:SYSSIZE
    diffOp(k,k) = linop.diff(DOMAIN, DIFFORDER(k));
end
t0 = tt(1);

% The vertical scale of the intial condition:
vscl = get(u0, 'vscale');

% Initial condition:
uCurrent = u0;
% storage
if ( SYSSIZE == 1 )
    uOut(1) = uCurrent;
else
    uOut{1} = uCurrent;
end

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
        currentLength = length(simplify(uCurrent, tol));
        pref.techPrefs.minPoints = max(2*currentLength, 9);
        chebfun( @(x) vscl + oneStep(x), DOMAIN, pref);
    else
        % Non-adaptive in space:
        currentLength = optN;
        oneStep(chebpts(optN, DOMAIN));
    end
    
    % Get chebfun of solution from this time chunk:
    uCurrent = chebfun(unew, DOMAIN);
    
    % Store in uOut:
    if ( SYSSIZE == 1 )
        % TODO?
        out = [uOut uCurrent];
%         uOut(nt+1) = uCurrent;
    else
        for k = 1:SYSSIZE
            uOut{nt+1} = uCurrent;
        end
    end
    
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
clear global GLOBX
clear global DMAT
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
            
            % Set the global variable x
            GLOBX = x;      
            
            % Compute the new differentiation matrices:
            makeDMAT(n);

            % Linear constraints:
            B = matrix([bc.left.op ; bc.right.op], n);
            
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
            F = pdeFun(U, t, x);
            F = M*F(:);
            
            % Get the algebraic right-hand sides: (may be time-dependent)
            for l = 1:numel(rhs)
                if ( isa(rhs{l}, 'function_handle') )
                    q(l, 1) = feval(rhs{l}, t);
                else
                    q(l, 1) = rhs{l};
                end
            end
            
            % Replacements for the BC algebraic conditions:
            F(rows) = B*U(:) - q;
            
            % Replacements for the nonlinear BC conditions:
            indx = 1:length(nllbc);
            if ( funFlagL )
                tmp = feval(nlbcsl, U, t, x);
                if ( size(tmp, 1) ~= n )
                    tmp = reshape(tmp, n, numel(tmp)/n);
                end
                F(rows(indx)) = tmp(1, :);
            else
                counter = 0;
                for l = 1:length(nllbc)
                    counter = counter + 1;
                    tmp = feval(nlbcs{counter}, U, t, x);
                    F(rows(l)) = tmp(1)-q(l);
                end
            end
            indx = numel(rhs) + 1 - nlrbc;
            if ( funFlagR )
                tmp = feval(nlbcsr, U, t, x);
                if ( size(tmp, 1) ~= n )
                    tmp = reshape(tmp, n, numel(tmp)/n);
                end
                F(rows(indx)) = fliplr(tmp(end, :));
            else
                for l = numel(rhs) + 1 - nlrbc
                    counter = counter + 1;
                    tmp = feval(nlbcs{counter}, U, t, x);
                    F(rows(l)) = tmp(end)-q(l);
                end
            end
            
            % Reshape to back to a single column:
            F = F(:);
        end
    end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DIFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The differential operators
function up = Diff(u, k)
% Computes the k-th derivative of u using Chebyshev differentiation
% matrices defined by barymat.

global DIFFORDER DMAT

% Assume first-order derivative:
if ( nargin == 1 )
    k = 1; 
end

% For finding the diff order of the RHS operator:
if ( any(isnan(u)) )
    idx = find(isnan(u), 1);
    if ( isempty(DIFFORDER) )
        DIFFORDER(idx) = k;
    else
        DIFFORDER(idx) = max(DIFFORDER(idx), k); 
    end
    up = u;
    return
end

% Trivial scalar and zero cases:
if ( size(u, 1) == 1 || ~any(u) )
    up = 0*u;
    return
end

% Find the derivative by muliplying by the kth-order differentiation matrix
up = DMAT{k}*u;

end

function makeDMAT(n)
    global DIFFORDER DOMAIN DMAT
    % Make the diffmats for this n:
    c = 2/diff(DOMAIN([1 end])); % Interval scaling
    DMAT = cell(max(DIFFORDER), 1);
    for kk = max(DIFFORDER):-1:1
        % TODO: Can this be done more efficiently?
        DMAT{kk} = c^kk*diffmat(n, kk);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SUM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The differential operators
function I = Sum(u, a, b)
% Computes the integral of u using clenshaw-curtis nodes and weights
% (which are stored for speed).

global GLOBX DOMAIN
persistent W
if ( isempty(W) )
    W = {};
end

% For finding the order of the RHS:
if ( any(isnan(u)) )
    I = u;
    return
end

N = length(u);

% Deal with the 3 args case. This can be integrating a sub-domain or
% indefinite integration. (Or integrating the whole domain...)
if ( nargin == 3 )
    x = GLOBX;
    if ( length(b) > 1 )
        if ( ~all(b == x) )
            error('CHEBFUN:pde15s:sumb', ...
                ['Limits in sum must be scalars or the indep space var ', ...
                 '(typically ''x'').']);
        elseif ( a < x(1) )
            error('CHEBFUN:pde15s:sumint', ...
                'Limits of integration outside of domain.');
        end
        
        I = Cumsum(u);
        I = I - bary(a, I, x);
        return
    elseif ( length(a) > 1 )
        if ( ~all(a == x) )
            error('CHEBFUN:pde15s:suma', ...
                ['Limits in sum must be scalars or the indep space var ', ...
                 '(typically ''x'').']);
        elseif ( b > x(end) )
            error('CHEBFUN:pde15s:sumint', ...
                'Limits of integration outside of domain.');
        end
        I = Cumsum(u);
        I = bary(b, I, x) - I;
        return
    elseif ( a ~= x(1) || b ~= x(end) )
        if ( a < x(1) || b > x(end) )
            error('CHEBFUN:pde15s:sumint', ...
            'Limits of integration outside of domain.');
        end
        I = Cumsum(u);
        I = bary(b, I, x) - bary(a, I, x);
        return
    else
        % Sum(u, a, b) is the same as below!
    end
end

% Retrieve or compute weights::
if ( N > 5 && numel(W) >= N && ~isempty(W{N}) )
    % Weights are already in storage!
else
    c = diff(DOMAIN)/2; % Interval scaling.
    W{N} = c*chebtech2.quadwts(N);
end

% Find the sum by muliplying by the weights vector:
I = W{N}*u;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CUMSUM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The differential operators
function U = Cumsum(u)
% Computes the indefinite integral of u.

global DOMAIN
persistent CMAT

% For finding the order of the RHS:
if ( any(isnan(u)) )
    U = u;
    return
end

N = length(u);
if ( N == 1 )
    U = u; 
    return
end

% Compute cumsum matrix:
if ( numel(CMAT) ~= N )
    c = diff(DOMAIN)/2; % Interval scaling.
    CMAT = cumsummat(N);
end

% Compute the indefinite integral:
U = CMAT*u;

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%  MISC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outfun = parseFun(inFun, flag)
% Rewrites the input function handle to call the right DIFF, SUM, methods, etc,
% and convert the input @(u, v, w, x, t) to @(U, x, t).

global SYSSIZE

Nin = nargin(inFun);
tmp = NaN(1, SYSSIZE);
% Determine number of operators, (i.e. diff, sum, cumsum) present in infun:
k = 1; Nops = [];
opsList = {@Diff, @Sum, @Cumsum};
while ( k < 4 && isempty(Nops) )
    tmp2 = repmat({tmp}, 1, nargin(inFun)-k);
    ops = opsList(1:k);
    try
        inFun(tmp2{:}, ops{:});
        Nops = k;
    catch ME
        try 
            % Try setting t to a scalar.
            tmp2{SYSSIZE+1} = NaN;
            inFun(tmp2{:}, ops{:});
            Nops = k;
        catch ME2
            %
        end
        % 
    end
    k = k+1;
end
if ( isempty(Nops) )
    error('CHEBFUN:pde15s:inputs', 'Unable to parse input function.');
end

% Convert inFun to accept quasimatrix inputs and remove ops from fun_handle:
ops = opsList(1:Nops);

% We don't accept only time or space as input args (i.e., both or nothing).
Nind = Nin - Nops - SYSSIZE;
if ( Nind == 0 )
    outfun = @(u, t, x) conv2cell(inFun, u, ops{:});
elseif ( Nind == 2 )
    outfun = @(u, t, x) conv2cell(inFun, u, t, x, ops{:});
else
    error('CHEBFUN:pde15s:inputs_ind', ['Incorrect number of independant ' ...
        'variables in input function. (Must be 0 or 2).']);
end

    function newFun = conv2cell(oldFun, u, varargin)
        % This function allows the use of different variables in the anonymous
        % function, rather than using the clunky quasi-matrix notation.
        
        % Note. This is faster than mat2cell!
        tmpCell = cell(1, SYSSIZE);
        if ( isa(u, 'chebfun') )
            for qk = 1:SYSSIZE
                tmpCell{qk} = u(qk);
            end
        else
            for qk = 1:SYSSIZE
                tmpCell{qk} = u(:, qk);
            end
        end
        newFun = oldFun(tmpCell{:}, varargin{:});
    end

end

function getDIFFORDER(pdeFun)
% Extract the DIFFORDER by evaluating the operator at some NaN values.
global DIFFORDER SYSSIZE
% Find the order of each of the variables:
DIFFORDER = zeros(1,SYSSIZE);
% tmp = repmat({zeros(1, SYSSIZE).'}, 1, SYSSIZE)
tmp = zeros(SYSSIZE);
for k = 1:SYSSIZE
    tmp(k,k) = NaN;
    pdeFun(tmp, 0, 0);
    tmp(k,k) = 0;
end
end

function out = dealWithStructInput(in)
if ( isstruct(in) )
    if ( numel(in) == 1)
        warning('PDE15S no longer supports struct inputs for bc.left and bc.right.')
        if ( isfield(in, 'val') )
            out = {in.op, in.val};
        else
            out = in.op;
        end
    else
        error('PDE15S no longer supports struct inputs for bc.left and bc.right.')
    end
else
    out = in;
end
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

function Q = cumsummat(N)

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
% CUMSUMMAT  Chebyshev integration matrix.
% Q = CUMSUMMAT(N) is the matrix that maps function values at N Chebyshev
% points to values of the integral of the interpolating polynomial at
% those points, with the convention that the first value is zero.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

N = N-1;

persistent cache    % stores computed values for fast return
if isempty(cache), cache = {}; end    % first call

if length(cache) >= N && ~isempty(cache{N})
    Q = cache{N};
    return
else
    cache{N} = [];
end

% Matrix mapping coeffs -> values.
T = cp2cdm(N);

% Matrix mapping values -> coeffs.
Tinv = cd2cpm(N);

% Matrix mapping coeffs -> integral coeffs. Note that the highest order
% term is truncated.
k = 1:N;
k2 = 2*(k-1);  k2(1) = 1;  % avoid divide by zero
B = diag(1./(2*k),-1) - diag(1./k2,1);
v = ones(N,1); v(2:2:end) = -1;
B(1,:) = sum( diag(v)*B(2:N+1,:), 1 );
B(:,1) = 2*B(:,1);

Q = T*B*Tinv;
Q(1,:) = 0;  % make exact
cache{N} = Q;

end

function T = cp2cdm(N)
% Values of Cheb. polys at Cheb nodes, x(n)=-cos(pi*n/N).
theta = pi*(N:-1:0)'/N;
T = cos( theta*(0:N) );
end

function C = cd2cpm(N)
% Three steps: Double the data around the circle, apply the DFT matrix,
% and then take half the result with 0.5 factor at the ends.
theta = (pi/N)*(0:2*N-1)';
F = exp( -1i*theta*(0:2*N-1) );  % DFT matrix
rows = 1:N+1;  % output upper half only
% Impose symmetries on data and coeffs.
C = real( [ F(rows,N+1) F(rows,N:-1:2)+F(rows,N+2:2*N) F(rows,1) ] );
C = C/N;  C([1 N+1],:) = 0.5*C([1 N+1],:);
end
