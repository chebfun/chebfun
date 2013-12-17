function varargout = pde15s( pdeFun, tt, u0, bc, varargin)
% PDE15S  Solve PDEs using the chebfun system.
%
% UU = PDE15s(PDEFUN, TT, U0, BC) where PDEFUN is a handle to a function
% with arguments u, t, x, and D, TT is a vector, U0 is a chebfun, and BC is
% a chebop boundary condition structure will solve the PDE
% dUdt = PDEFUN(UU, t, x) with the initial condition U0 and boundary
% conditions BC over the time interval TT.
%
% PDEFUN should take the form @(U1, U2, ..., UN, T, X, D, S, C), where U1, ..., UN
% are the unknown dependent variables to be solved for, T is time, X is
% space, D is the differential operator, S is the definite integral
% operator (i.e., 'sum'), and C the indefinite integral operator (i.e.,
% 'cumsum').
%
% For equations of one variable, UU is output as a quasimatrix, where UU(:, k)
% is the solution at TT(k). For systems, the solution is returned as a
% cell array of quasimatrices.
%
% Example 1: Nonuniform advection
%   x = chebfun('x', [-1 1]);
%   u = exp(3*sin(pi*x));
%   f = @(u, t, x, diff) -(1+0.6*sin(pi*x)).*diff(u);
%   uu = pde15s(f, 0:.05:3, u, 'periodic');
%   surf(u, 0:.05:3)
%
% Example 2: Kuramoto-Sivashinsky
%   d = domain(-1, 1);
%   x = chebfun('x');
%   I = eye(d); D = diff(d);
%   u = 1 + 0.5*exp(-40*x.^2);
%   bc.left = struct('op', {I, D}, 'val', {1, 2});
%   bc.right = struct('op', {I, D}, 'val', {1, 2});
%   f = @(u, diff) u.*diff(u)-diff(u, 2)-0.006*diff(u, 4);
%   uu = pde15s(f, 0:.01:.5, u, bc);
%   surf(u, 0:.01:.5)
%
% Example 3: Chemical reaction (system)
%    x = chebfun('x', [-1 1]);
%    u = [ 1-erf(10*(x+0.7)) , 1 + erf(10*(x-0.7)) , 0 ];
%    f = @(u, v, w, diff)  [ .1*diff(u, 2) - 100*u.*v , ...
%                         .2*diff(v, 2) - 100*u.*v , ...
%                         .001*diff(w, 2) + 2*100*u.*v ];
%    bc = 'neumann';
%    uu = pde15s(f, 0:.1:3, u, bc);
%    mesh(uu{3})
%
% See chebfun/examples/pde15s_demos.m and chebfun/examples/pde_systems.m
% for more examples.
%
% UU = PDE15s(PDEFUN, TT, U0, BC, OPTS) will use nondefault options as
% defined by the structure returned from OPTS = PDESET.
%
% UU = PDE15s(PDEFUN, TT, U0, BC, OPTS, N) will not adapt the grid size
% in space. Alternatively OPTS.N can be set to the desired size.
%
% [TT UU] = PDE15s(...) returns also the time chunks TT.
%
% There is some support for nonlinear and time-de[pendent boundary
% conditions, such as
%    BC.LEFT = @(u, t, x, diff) diff(u) + u.^2 - (1+2*sin(10*t));
%    BC.RIGHT = struct( 'op', 'dirichlet', 'val', @(t) .1*sin(t));
% with the input format being the same as PDEFUN described above.
%
% See also PDESET, ODE15S, CHEBOP/PDE15S.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

global ORDER GLOBX DMat
ORDER = 0; % Initialise to zero
GLOBX = [];
DMat = {};

% if ( nargin < 4 )
%     error('CHEBFUN:pde15s:argin', 'pde15s requires a minimum of 4 inputs.');
% end

% Default options:
tol = 1e-6;             % 'eps' in chebfun terminology
doPlot = 1;             % plot after every time chunk?
doHold = 0;             % hold plot?
plotOpts = '-';         % Plot Style

% Parse the variable inputs
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

% PDE solver options
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

% Experimental feature for coupled ode/pde systems
if ( isfield(opt, 'PDEflag') )
    pdeFlag = opt.PDEflag;
else
    pdeFlag = true;
end

% % Determine which figure to plot to (for CHEBGUI) and set default display values
% % for variables.
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
    varnames = opt.handles.varnames;
    xLabel = opt.handles.indVarName{1};
    tlabel = opt.handles.indVarName{2};
else
    varnames = 'u';
    xLabel = 'x';
    tlabel = 't';
end

% Parse plotting options
indx = strfind(plotOpts, ', ');
tmpOpts = cell(numel(indx)+1, 1);
k = 0; j = 1;
while ( k < numel(plotOpts) )
    k = k+1;
    sk = plotOpts(k);
    if strcmp(sk, ', ')
        tmpOpts{j} = plotOpts(1:k-1);
        plotOpts(1:k) = [];
        j = j+1;
        k = 0;
    end
end
tmpOpts{j} = plotOpts;
plotOpts = tmpOpts;
for k = 1:numel(plotOpts)
    if strcmpi(plotOpts{k}, 'linewidth') || strcmpi(plotOpts{k}, 'MarkerSize')
        plotOpts{k+1} = str2double(plotOpts{k+1});
    end
end

% ODE tolerances:
% (AbsTol and RelTol must be <= Tol/10)
aTol = odeget(opt, 'AbsTol', tol/10);
rTol = odeget(opt, 'RelTol', tol/10);
if ( isnan(optN) )
    aTol = min(aTol, tol/10);
    rTol = min(rTol, tol/10);
end
opt.AbsTol = aTol;
opt.RelTol = rTol;

% Check for (and try to remove) piecewise initial conditions
u0Trans = u0(1).isTransposed;
% Get u0trans as a cell if u0 is a quasimatrx
if ( u0Trans )
    u0 = transpose(u0);
end

for k = 1:numel(u0)
    if ( numel(u0(k).funs) > 1 )
        u0(k) = merge(u0(k), 'all', 1025, tol);
        if ( u0(k).nfuns > 1 )
            error('CHEBFUN:pde15s:piecewise', ...
                'Piecewise initial conditions are not supported');
        end
    end
end
% Simplify initial condition to tolerance or fixed size in optN
if ( isnan(optN) )
    u0 = simplify(u0);
else
    for k = 1:numel(u0)
        u0(k).funs{1}.onefun = prolong(u0(k).funs{1}.onefun, optN);
    end
end

% Get the domain and the independent variable 'x':
d = domain(u0);
xd = chebfun(@(x) x, d);

% These are used often.
Z = linop.zeros(d);
I = linop.eye(d);
D = linop.diff(d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% parse inputs to pdefun %%%%%%%%%%%%%%%%%%%%%%%%%%%%
sysSize = min(size(u0));            % Determine the size of the system
pdeFun = parsefun(pdeFun, sysSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% parse boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%
% Some error checking on the bcs
if ( ischar(bc) && (strcmpi(bc, 'neumann') || strcmpi(bc, 'dirichlet')) )
    bc = struct( 'left', bc, 'right', bc);
elseif ( iscell(bc) && numel(bc) == 2 )
    bc = struct( 'left', bc{1}, 'right', bc{2});
end

% Shorthand bcs - all neumann or all dirichlet
if ( isfield(bc, 'left') && (ischar(bc.left) || (iscell(bc.left) && ischar(bc.left{1}))) )
    if ( iscell(bc.left) )
        v = bc.left{2};
        bc.left = bc.left{1};
    else
        v = 0;
    end
    if ( strcmpi(bc.left, 'dirichlet') )
        A = I;
    elseif ( strcmpi(bc.left, 'neumann') )
        A = D;
    end
    op = cell(1, sysSize);
    
    for k = 1:sysSize
        op{k} = [repmat(Z, 1, k-1) A repmat(Z, 1, sysSize-k)];
    end
    bc.left = struct('op', op, 'val', repmat({v}, 1, sysSize));
end
if ( isfield(bc, 'right') && (ischar(bc.right) || (iscell(bc.right) && ischar(bc.right{1}))) )
    if ( iscell(bc.right) )
        v = bc.right{2};
        bc.right = bc.right{1};
    else
        v = 0;
    end
    if ( strcmpi(bc.right, 'dirichlet') )
        A = I;
    elseif ( strcmpi(bc.right, 'neumann') )
        A = D;
    end
    op = cell(1, sysSize);
    for k = 1:sysSize
        op{k} = [repmat(Z, 1, k-1) A repmat(Z, 1, sysSize-k)];
    end
    bc.right = struct('op', op, 'val', repmat({v}, 1, sysSize));
end

if ( isfield(bc, 'left') && ~isfield(bc, 'right') )
    bc.right = [];
elseif ( isfield(bc, 'right') && ~isfield(bc, 'left') )
    bc.left = [];
end

% Sort out left boundary conditions
nllbc = []; nlbcs = {}; GLOBX = 1; funFlagL = false; rhs = {};
% 1) Deal with the case where bc is a function handle vector
if ( isfield(bc, 'left') && numel(bc.left) == 1 && isa(bc.left, 'function_handle') )
    op = parsefun(bc.left, sysSize);
    sop = size(op(ones(1, sysSize), 0, mean(d.ends)));
    nllbc = 1:max(sop);
    bc.left = struct( 'op', [], 'val', []);
    % Dummy entries (Worked out naively. AD information may be used below).
    for k = nllbc
        if sysSize == 1
            bc.left(k).op = repmat(I, 1, sysSize);
        else
            bc.left(k).op = [repmat(Z, 1, k-1) I repmat(Z, 1, sysSize-k)];
        end
        bc.left(k).val = 0;
    end
    rhs = num2cell(zeros(1, max(sop)));
    nlbcsl = op;
    funFlagL = true;
    % 2) Deal with other forms of input
elseif ( isfield(bc, 'left') && numel(bc.left) > 0 )
    if ( isa(bc.left, 'linop') || iscell(bc.left) )
        bc.left = struct( 'op', bc.left);
    elseif ( isnumeric(bc.left) )
        bc.left = struct( 'op', I, 'val', bc.left);
    end
    % Extract nonlinear conditions
    rhs = cell(numel(bc.left), 1);
    for k = 1:numel(bc.left)
        opk = bc.left(k).op;
        rhs{k} = 0;
        
        % Numerical values
        if ( isnumeric(opk) && sysSize == 1 )
            bc.left(k).op = repmat(I, 1, sysSize);
            bc.left(k).val = opk;
        end
        
        % Function handles
        if ( isa(opk, 'function_handle') )
            nllbc = [nllbc k];             % Store positions
            nlbcs = [nlbcs {parsefun(opk)}];
            % Dummy entries (Worked out naively. AD information may be used below).
            bc.left(k).op = [repmat(Z, 1, k-1) I repmat(Z, 1, sysSize-k)];
        end
        
        % Remove 'vals' from bc and construct cell of rhs entries
        if ( isfield(bc.left(k), 'val') && ~isempty(bc.left(k).val) )
            rhs{k} = bc.left(k).val;
        end
        bc.left(k).val = 0;  % remove function handles
    end
end

% Sort out right boundary conditions
nlrbc = []; numlbc = numel(rhs); funFlagR = false;
% 1) Deal with the case where bc is a function handle vector
if ( isfield(bc, 'right') && numel(bc.right) == 1 && isa(bc.right, 'function_handle') )
    op = parsefun(bc.right, sysSize);
    sop = size(op(ones(1, sysSize), 0, mean(d.ends)));
    nlrbc = 1:max(sop);
    bc.right = struct( 'op', [], 'val', []);
    % Dummy entries (Worked out naively. AD information may be used below).
    for k = nlrbc
        if ( sysSize == 1 )
            bc.right(k).op = repmat(I, 1, sysSize);
        else
            bc.right(k).op = [repmat(Z, 1, k-1) I repmat(Z, 1, sysSize-k)];
        end
        bc.right(k).val = 0;
    end
    rhs = [rhs num2cell(zeros(1, max(sop)))];
    nlbcs = op;
    funFlagR = true;
    % 2) Deal with other forms of input
elseif ( isfield(bc, 'right') && numel(bc.right) > 0 )
    if ( isa(bc.right, 'linop') || isa(bc.right, 'cell') )
        bc.right = struct( 'op', bc.right, 'val', 0);
    elseif ( isnumeric(bc.right) )
        bc.right = struct( 'op', I, 'val', bc.right);
    end
    for k = 1:numel(bc.right)
        opk = bc.right(k).op;
        rhs{numlbc+k} = 0;
        if isnumeric(opk) && sysSize == 1
            bc.right(k).op = I;
            bc.right(k).val = opk;
        end
        if ( isa(opk, 'function_handle') )
            nlrbc = [nlrbc k];
            nlbcs = [nlbcs {parsefun(opk)}];
            % Dummy entries (Worked out naively. AD information may be used below).
            bc.right(k).op = [repmat(Z, 1, k-1) I repmat(Z, 1, sysSize-k)];
        end
        if ( isfield(bc.right(k), 'val') && ~isempty(bc.right(k).val) )
            rhs{numlbc+k} = bc.right(k).val;
        end
        bc.right(k).val = 0;
    end
end

t0 = tt(1);

%%%%%%%%%% Support for user-defined mass matrices and coupled BVP-PDEs! %%%%%%%
if ( ~isempty(opt.Mass) || ~all(pdeFlag) )
    userMassSet = true;
    userMassMatrix = [];
    if ( ~all(pdeFlag) )
        for k = 1:numel(pdeFlag)
            if ( pdeFlag(k) )
                A = I;
            else
                A = Z;
            end
            userMassMatrix = [userMassMatrix ; repmat(Z, 1, k-1) A repmat(Z, 1, sysSize-k)];
        end
    end
else
    userMassSet = false;
end

% This is needed inside the nested function onestep()
for k = 1:numel(sysSize)
    diffOp{k} = linop.diff(d, ORDER(k));
end
diffOp = [diffOp{:}];

% The vertical scale of the intial condition
vscl = get(u0, 'vscale');

% Plotting setup
if ( doPlot )
    if ( ~guiFlag )
        cla, shg
    end
    set(gcf, 'doublebuf', 'on');
    if ( isempty(get(u0(:, 1), 'impulses')) )
        for k = 1:numel(u0);
            u0(:, k) = set(u0(:, k), 'impulses', get(u0(:, k), 'domain'));
        end
    end
    plot(u0, plotOpts{:});
    if ( doHold )
        ish = ishold;
        hold on
    end
    if ( ~isempty(YLim) )
        ylim(YLim);
    end
    % Axis labels
    xlabel(xLabel);
    if ( numel(varnames) > 1 )
        legend(varnames);
    else
        ylabel(varnames);
    end
    % Determines whether grid is on
    if ( gridOn )
        grid on
    end
    drawnow
end

% initial condition
uCurrent = u0;
% storage
if ( sysSize == 1 )
    uOut(1) = uCurrent;
else
    uOut{1} = uCurrent;
end

% initialise variables for onestep()
B = []; q = []; rows = []; M = []; n = [];

% Set the preferences
pref = chebpref;
pref.techPrefs.eps = tol;
pref.refinementFunction = 'resampling'; 
pref.enableBreakpointDetection = 0;
pref.techPrefs.sampleTest = 0; 
pref.enableSingularityDetection = 0; 

% Begin time chunks
for nt = 1:length(tt)-1
    
    % Solve one chunk:
    if ( isnan(optN) )
        % Size of current length:
        currentLength = length(simplify(uCurrent, tol));
        pref.techPrefs.minPoints = max(currentLength,9);
        chebfun( @(x) vscl + onestep(x), d, pref);
    else
        % Non-adaptive in space:
        onestep(chebpts(optN, d));
    end
    
    % Get chebfun of solution from this time chunk:
    uCurrent = chebfun(unew, d);
    
    % Store in uOut:
    if ( sysSize == 1 )
        uOut(nt+1) = uCurrent;
    else
        for k = 1:sysSize
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
        if ( numel(varnames) > 1 )
            legend(varnames);
        else
            ylabel(varnames);
        end
        % Determines whether grid is on
        if ( gridOn )
            grid on
        end
        title(sprintf('%s = %.3f,  len = %i', tlabel, tt(nt+1), currentLength)), drawnow
    elseif ( guiFlag )
        drawnow
    end
    
    if ( guiFlag )
        % Interupt comutation if stop or pause  button is pressed in the GUI.
        if ( strcmp(get(solveButton, 'String'), 'Solve') )
            tt = tt(1:nt+1);
            if sysSize == 1,
                uOut = uOut(1:nt+1);
            else
                for k = 1:sysSize
                    uOut{k} = uOut{k}(1:nt+1);
                end
            end
            break
        elseif ( strcmp(get(clearButton, 'String'), 'Continue') )
            defaultlinewidth = 2;
            axes(axesNorm)
            if ( ~iscell(uOut) )
                waterfall(uOut(1:nt+1), tt(1:nt+1), 'simple', 'linewidth', defaultlinewidth)
                xlabel(xLabel), ylabel(tlabel), zlabel(varnames)
            else
                cols = get(0, 'DefaultAxesColorOrder');
                for k = 1:numel(uOut)
                    plot(0, NaN, 'linewidth', defaultlinewidth, 'color', cols(k, :)), hold on
                end
                legend(varnames);
                for k = 1:numel(uOut)
                    waterfall(uOut{k}, tt(1:nt+1), 'simple', 'linewidth', defaultlinewidth, 'edgecolor', cols(k, :)), hold on
                    xlabel(xLabel), ylabel(tlabel)
                end
                view([322.5 30]), box off, grid on, hold off
            end
            axes(axesSol)
            waitfor(clearButton, 'String');
        end
    end
end

if ( doPlot && doHold && ~ish )
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
        error('CHEBFUN:pde15s:output', ...
            'pde15s may only have a maximum of two outputs.');
end

clear global ORDER
clear global GLOBX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    ONESTEP   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs the result of one time chunk at fixed discretization
    function U = onestep(x)
        
        if ( length(x) == 2 )
            U = [0;0]; 
            return
        end
        
        % Evaluate the chebfun at discrete points
        U0 = feval(uCurrent, x);
        % This depends only on the size of n. If this is the same, reuse
        if ( isempty(n) || n ~= length(x) )
            
            % The new discretisation length
            n = length(x);
            
            % Set the global variable x
            GLOBX = x;      
            
            % Make the diffmats for this n:
            c = 2/diff(x([1 end])); % Interval scaling
            DMat = cell(max(ORDER), 1);
            for m = max(ORDER):-1:1
                DMat{m} = c^m*diffmat(n, m);
            end
            
            % See what the boundary replacement actions will be.
%             [ignored, B, q, rows] = matrix( addbc(diffOp, bc), n, 'oldschool');
            % TODO!
%             B = zeros(1, n); B([1, end]) = [-1, 1];
%             q = 0;
%             rows = 1;
            
            In = eye(n);
            Dn = diffmat(n);
            Zn = zeros(2,n);
%             BD = [In([1 end],:) Zn Zn ; Zn In([1 end],:) Zn ; Zn Zn In([1 end],:)];
            BN = [Dn([1 end],:) Zn Zn ; Zn Dn([1 end],:) Zn ; Zn Zn Dn([1 end],:)];
            B = [BN];
            
            q = zeros(6, 1);
            rows = [];
            for k = 1:sysSize
                rows = [rows, (k-1)*n + [1 n]];
            end
            
            % Mass matrix is I except for algebraic rows for the BCs.
            M = speye(sysSize*n);    
            M(rows, :) = 0;
            
            % Multiply by user-defined mass matrix
            if ( userMassSet ) 
                M = feval(userMassMatrix, n)*M; 
            end
            
        end
        
        % ODE options: (mass matrix)
        opt2 = odeset(opt, 'Mass', M, 'MassSingular', 'yes', ...
            'InitialSlope', odeFun(tt(nt), U0), 'MStateDependence', 'none');
        
        % Solve ODE over time chunk with ode15s:
        [ignored, U] = ode15s(@odeFun, tt(nt:nt+1), U0, opt2);
        
        % Reshape solution:
        U = reshape(U(end, :).', n, sysSize);
        
        % The solution we'll take out and store:
        unew = U;
        
        % Collapse systems to single chebfun for constructor (is addition right?)
        U = sum(U, 2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    ODEFUN   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % This is what ode15s calls.
        function F = odeFun(t, U)
            
            % Reshape to n by syssize
            U = reshape(U, n, sysSize);
            
            % Evaluate the PDEFUN
            F = pdeFun(U, t, x);
            
            % Get the algebraic right-hand sides (may be time-dependent)
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
                j = 0;
                for kk = 1:length(nllbc)
                    j = j + 1;
                    tmp = feval(nlbcs{j}, U, t, x);
                    F(rows(kk)) = tmp(1)-q(kk);
                end
            end
            indx = numel(rhs)+1-nlrbc;
            if ( funFlagR )
                tmp = feval(nlbcsr, U, t, x);
                if ( size(tmp, 1) ~= n )
                    tmp = reshape(tmp, n, numel(tmp)/n);
                end
                F(rows(indx)) = fliplr(tmp(end, :));
            else
                for kk = numel(rhs)+1-nlrbc
                    j = j + 1;
                    tmp = feval(nlbcs{j}, U, t, x);
                    F(rows(kk)) = tmp(end)-q(kk);
                end
            end
            
            % Reshape to single column
            F = F(:);
            
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   DIFF   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The differential operators
function up = Diff(u, k)
% Computes the k-th derivative of u using Chebyshev differentiation
% matrices defined by barymat.

global ORDER DMat

% Assume first-order derivative
if ( nargin == 1 )
    k = 1; 
end

if ( isa(u, 'chebfun') )
    up = diff(u, k); 
    return
end

% For finding the order of the RHS
if ( any(isnan(u)) )
    idx = find(isnan(u), 1);
    if ( isempty(ORDER) )
        ORDER(idx) = k;
    else
        ORDER(idx) = max(ORDER(idx), k); 
    end
    up = u;
    return
end

if ( size(u, 1) == 1 )
    up = 0*u;
    return
end

% Find the derivative by muliplying by the kth-order differentiation matrix
up = DMat{k}*u;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   SUM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The differential operators
function I = Sum(u, a, b)
% Computes the integral of u using clenshaw-curtis nodes and weights
% (which are stored for speed).

global GLOBX
persistent W
if ( isempty(W) )
    W = {};
end

if ( isa(u, 'chebfun') )
    if ( nargin == 1 )
        I = sum(u);
    else
        I = sum(u, a, b);
    end
    return
end

% For finding the order of the RHS
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
                'Limits in sum must be scalars or the indep space var (typically ''x'').');
        elseif ( a < x(1) )
            error('CHEBFUN:pde15s:sumint', 'Limits of integration outside of domain.');
        end
        
        I = Cumsum(u);
        I = I - bary(a, I, x);
        return
    elseif ( length(a) > 1 )
        if ( ~all(a == x) )
            error('CHEBFUN:pde15s:suma', ...
                'Limits in sum must be scalars or the indep space var (typically ''x'').');
        elseif ( b > x(end) )
            error('CHEBFUN:pde15s:sumint', 'Limits of integration outside of domain.');
        end
        I = Cumsum(u);
        I = bary(b, I, x) - I;
        return
    elseif ( a ~= x(1) || b ~= x(end) )
        if ( a < x(1) || b > x(end) )
            error('CHEBFUN:pde15s:sumint', 'Limits of integration outside of domain.');
        end
        I = Cumsum(u);
        I = bary(b, I, x) - bary(a, I, x);
        return
    end
end

% Retrieve or compute weights:
if ( N > 5 && numel(W) >= N && ~isempty(W{N}) )
    % Weights are already in storage!
else
    x = GLOBX;
    c = diff(x([1 end]))/2;
    W{N} = c*weights2(N);
end

% Find the sum by muliplying by the weights vector:
I = W{N}*u;

end

function w = weights2(n) % 2nd kind
% Jörg Waldvogel, "Fast construction of the Fejér and Clenshaw-Curtis
% quadrature rules", BIT Numerical Mathematics 43 (1), p. 001-018 (2004).
if ( n == 1 )
    w = 2;
else
    m = n-1;
    c = zeros(1, n);
    c(1:2:n) = 2./[1 1-(2:2:m).^2 ];
    f = real(ifft([c(1:n) c(m:-1:2)]));
    w = [f(1) 2*f(2:m) f(n)];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   CUMSUM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The differential operators
function U = Cumsum(u)
% Computes the indefinite integral of u.

global GLOBX

if ( isa(u, 'chebfun') )
    U = cumsum(u); 
    return
end

% For finding the order of the RHS
if ( any(isnan(u)) )
    U = u;
    return
end

N = length(u);
if ( N == 1 )
    U = u; 
    return
end

% Compute matrix.
x = GLOBX;
c = diff(x([1 end]))/2;
C = cumsummat(N);

% Find the indefinite integral by muliplying cumsum matrix:
U = c*(C*u);
end

function outfun = parsefun(inFun, sysSize)

global ORDER

Nin = nargin(inFun);
tmp = NaN(1, sysSize);
% Number of operators, (i.e. diff, sum, cumsum) present in infun
% Also computes QUASIN through global variable in Diff.
k = 1; Nops = [];
opsList = {@Diff, @Sum, @Cumsum};
while ( k < 4 && isempty(Nops) )
    tmp2 = repmat({tmp}, 1, nargin(inFun)-(k+1));
    ops = opsList(1:k);
    inFun(tmp, tmp2{:}, ops{:});
    Nops = k;
    k = k+1;
end

% Find the order of each of the variables:
ORDER = zeros(1,sysSize);
tmp = repmat({zeros(1, sysSize)}, 1, nargin(inFun)-Nops);
for k = 1:sysSize
    tmp{k}(k) = NaN;
    inFun(tmp{:}, ops{:});
    tmp{k}(k) = 0;
end

if ( isempty(Nops) )
    error('CHEBFUN:pde15s:inputs', 'Unable to parse input function.');
end

Nind = Nin - Nops - sysSize;

% We don't accept only time or space as input args (both or nothing).
if ~(Nind == 0 || Nind == 2)
    error('CHEBFUN:pde15s:inputs_ind', ['Incorrect number of independant variables' ...
        ' in input function. (Must be 0 or 2).']);
end

% Convert inFun to accept quasimatrix inputs and remove ops from fun handle
ops = opsList(1:Nops);

if ( Nind == 0 )
    outfun = @(u, t, x) conv2cell(inFun, u, ops{:});
elseif ( Nind == 2 )
    outfun = @(u, t, x) conv2cell(inFun, u, t, x, ops{:});
end

    function newFun = conv2cell(oldFun, u, varargin)
        % This function allows the use of different variables in the anonymous
        % function, rather than using the clunky quasi-matrix notation.
        tmpCell = cell(1, sysSize);
        if ( isa(u, 'chebfun') )
            for qk = 1:sysSize
                tmpCell{qk} = u(qk);
            end
        else
            for qk = 1:sysSize
                tmpCell{qk} = u(:, qk);
            end
        end
        newFun = oldFun(tmpCell{:}, varargin{:});
    end

end
