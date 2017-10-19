function varargout = solveivp(N, rhs, pref, varargin)
%SOLVEIVP    Solve an IVP by reforming it to a first order system.
%
%   U = SOLVEIVP(N, RHS), where N is a CHEBOP and RHS is a CHEBMATRIX, CHEBFUN
%   or a vector of doubles attempts to solve the IVP
%
%       N(U) = RHS + boundary conditions specified by N
%
%   Observe that U = SOLVEIVP(N, RHS), where N specifies an initial/final-value
%   problem (IVP/FVP), has the same effect as U = N\RHS, but this method allows
%   greater flexibility than CHEBOP backslash, as described below. Problems are
%   determined to be an IVP/FVP as follows:
%       * N.LBC is non-empty, N.RBC and N.BC are empty => IVP.
%       * N.RBC is non-empty, N.LBC and N.BC are empty => FVP.
%   Otherwise, problems are considered to be boundary-value problems, and
%   U=N\RHS will in general have the same effect as U = SOLVEBVP(N, RHS).
%
%   If successful, the solution returned, U, is a CHEBFUN if N specifies a
%   scalar problem, and a CHEBMATRIX if N specifies a coupled systems of
%   ordinary differential equations. See note below on how to call the method
%   with multiple outputs. This method solves both linear and nonlinear problems
%   be automatically converting them to a coupled first-order system, which can
%   then be solved using MATLAB's built in solvers.
%
%   U = SOLVEIVP(N, RHS, PREF) is the same as above, using the preferences
%   specified by the CHEBOPPREF variable PREF.
%
%   [U, INFO] = SOLVEIVP(N, RHS, PREF) is the same as above, but also returns
%   the MATLAB struct INFO, which contains useful information about the solution
%   process. The fields of INFO are as follows (more may be added in future
%   versions):
%       SOLVER: The MATLAB solver used when solving the problem.
%
%   [U, V, ...] = SOLVEBVP(N, ...), where N specifies a coupled system of ODEs,
%   returns CHEBFUNs U, V, ... for individual solution components, rather than a
%   CHEBMATRIX.
%
%
%   Note 1: CHEBOP allows the RHS of coupled system of ODEs to be a scalar,
%   e.g., one can both call
%       N = chebop(@(x, u, v) [diff(u) - v.^2 ; u - diff(v)], [0 3]);
%       N.lbc = @(u, v) [u - 1 ; v + 1];
%       uv = solveivp(N, 0);
%   and
%       uv = solveivp(N, [0; 0]);
%
%
%   Note 2: The solver tries to construct global CHEBFUNs to represent the
%   solutions of ODEs if possible (however, breakpoints in the domain and
%   coefficients do get respected). Turn on the global CHEBFUN splitting option
%   if you wish to obtain solutions with further breakpoints, e.g.
%       % Solve van der Pol equation without and with splitting
%       N = chebop(@(t,u) diff(u,2)-5.*(1-u.^2).*diff(u)+u, [0 20]);
%       N.lbc = @(u) [u-2; diff(u)];
%       uNoSplit = N\0
%       % Turn on splitting with max length of each piece equal to 300
%       chebfunpref.setDefaults('splitting', true)
%       chebfunpref.setDefaults({'splitPrefs','splitLength'}, 300)
%       uSplit = N\0
%       % Turn splitting back off
%       chebfunpref.setDefaults('splitting', false)
%       
%
% See also: CHEBOP, CHEBOP/MLDIVIDE, CHEBOPPREF, CHEBOP/SOLVEBVP,
% CHEBFUN/ODE113, CHEBFUN/ODE15S, CHEBFUN/ODE45, CHEBFUN/ODESOL, TREEVAR. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Developer note:
%   U = SOLVEIVP(N, RHS, PREF, DISPLAYINFO) allows passing in a function handle
%   to a displaying method that is called after solution has been computed.
%   This allows separating the displaying process for regular CHEBOP use and
%   CHEBGUI. See chebop/displayIVPinfo() and chebgui/displayIVPinfo() for more
%   details.

%% SETUP
% Check inputs
if ( nargin < 2 )
    % Default right-hand side.
    rhs = 0;
elseif (isa(rhs, 'chebfun') || isa(rhs, 'chebmatrix') )
    % Ensure that the RHS lives on the same domain as the CHEBOP N.
    assert( (rhs.domain(1) == N.domain(1)) && ...
            (rhs.domain(end) == N.domain(end)), ...
            'CHEBFUN:CHEBOP:solveivp:domainMismatch', ...
            ['The domains of the CHEBOP and the right-hand side of the ' ...
            'equation do not match.']);
end

if ( nargin < 3 )
    % Create a default CHEBOPPREF if we didn't get one passed
    pref = cheboppref;
end

if ( nargin < 4 )
    % Use the default displayInfo() method if none was passed.
    displayInfo = @chebop.displayIVPinfo;
else
    % Otherwise, the displayInfo method was the first entry of VARARGIN.
    displayInfo = varargin{1};
end

% Check for unbounded domains, which we (currently) don't support:
if ( ~all(isfinite(N.domain)) )
    error('CHEBFUN:CHEBOP:solveivp:infDom', ...
        'Solving IVPs on unbounded intervals is not supported.');
end

% What solver do we want to use for the IVP?
solver = pref.ivpSolver;

% If pref.ivpSolver is set to a global method, we really should be calling
% CHEBOP/SOLVEBVP():
if ( strcmp(solver, 'values') || strcmp(solver, 'coeffs') || ...
        isempty(strfind(func2str(solver), 'chebfun.ode')) )
    [varargout{1:nargout}] = solvebvp(N, rhs, pref, varargin{:});
    info.solver = 'Global method';
    return
end

% Find out how many variables N operates on:
nVars = numVars(N);

% Check if we're working with a CHEBMATRIX syntax, e.g.
%   @(x,u) [diff(u{1}) + u{2}; ...]
% as we'll need different feval calls if that's the case.
if ( nargin(N) == 2 && nVars > 1 )
    cellArg = 1;
else
    cellArg = 0;
end

%% Convert to a first-order system

% We call the conversion methods of the TREEVAR class, the call depends on
% whether we're dealing with a system or not. We don't support integral
% equations, so listen out for warnings that the TREEVAR class throws if it
% encounters unsupported methods.
try
    [anonFun, varIndex, problemDom, coeffs, diffOrders] = ...
        treeVar.toFirstOrder(N.op, rhs, N.domain, nVars, cellArg);
catch ME
    % Did we encounter an unsupported method? If so, try to solve it globally:
    if ( ~isempty(regexp(ME.identifier, 'CHEBFUN:TREEVAR:.+:notSupported', ...
            'once')) )
        [varargout{1:nargout}] = solvebvp(N, rhs, pref, varargin{:});
        return
    else
        % Otherwise, an unexpected error occured, rethrow it.
        rethrow(ME);
    end
end
    
% Join all breakpoints, which can either be specified by the CHEBOP, or arise
% from discontinuous coefficients in the problem.
dom = union(N.domain, problemDom);

%% Obtain information about the initial conditions.

% Create a zero CHEBFUN for evaluating N.LBC or N.RBC. This gives us the correct
% values we need for passing as initial/final conditions to the ODE solver.
cheb0 = chebfun(@(x) 0*x, dom);

% Evaluate N.LBC or N.RBC:
if ( ~isempty(N.lbc) )
    % Create enough copies to allow to evaluate the initial condition for
    % systems. We use the DIFFORDER variable from the first order conversion
    % above for giving us information about the number of variables in the
    % problem:
    cheb0 = repmat({cheb0}, length(diffOrders), 1);
    % Evaluate the initial condition, depending on whether we're dealing with
    % CHEBMATRIX syntax or not:
    if ( cellArg && nargin(N.lbc) == 1 )
       bcEvalFun = N.lbc(cheb0); 
    else
        bcEvalFun = N.lbc(cheb0{:});
    end
    % Store the point at which the initial condition is evaluated (left
    % endpoint):
    evalPoint = dom(1);
    % The domain of the differential equation:
    odeDom = dom;
    % Store that we're solving an IVP.
    isIVP = 1;
else
    % Create enough copies to allow to evaluate the initial condition for
    % systems:
    cheb0 = repmat({cheb0}, length(diffOrders), 1);
    % Evaluate the final condition, depending on whether we're dealing with
    % CHEBMATRIX syntax or not:
    if ( cellArg && nargin(N.rbc) == 1 )
       bcEvalFun = N.rbc(cheb0); 
    else
        bcEvalFun = N.rbc(cheb0{:});
    end
    % Store the point at which the final condition is evaluated (right
    % endpoint):
    evalPoint = dom(end);
    % Flip the time domain of the problem, so that it will be solved as a
    % final-value problem.
    odeDom = fliplr(dom);
    % Store that we're solving an FVP.
    isIVP = 0;
end

% Check that the coefficients multiplying the highest order terms are
% non-vanishing at the end point of the domain, as the MATLAB ODE methods can't
% deal with such problems:
vanishingCoeffs = false;
for coeffCounter = 1:length(coeffs)
    coeff = coeffs{coeffCounter};
    if ( isa(coeff, 'chebfun') )
        vanishingCoeffs = vanishingCoeffs || ( coeff(evalPoint) == 0 );
    else
        vanishingCoeffs = vanishingCoeffs || ( coeff == 0 );
    end
end
if ( vanishingCoeffs )
    error('CHEBFUN:CHEBOP:solveivp:zeroCoeffs', ...
        ['Time stepping methods do not support vanishing coefficients at ' ...
        'the endpoint of the domain. Please use global methods instead.']);
end

% Evaluating N.LBC or N.RBC above gives CHEBFUN results. Evaluate the CHEBFUNs
% at the appropriate endpoint to obtain the scalar value for the initial/final
% condition. We need to negate the results so that we can pass them as
% appropriate to the MATLAB solvers below
if ( isa(bcEvalFun, 'chebfun') )
    % Scalar case
    initVals = -bcEvalFun(evalPoint);
else
    % Coupled systems case
    initVals = -cellfun(@feval, bcEvalFun.blocks, ...
        repmat({evalPoint}, length(bcEvalFun.blocks), 1) );
end

% The number of initial conditions passed should match the total diffOrders that
% appear in the problem:
assert(sum(diffOrders) == length(initVals), ...
    'CHEBFUN:CHEBOP:solveivp:numConditions', ['The number of initial/final ' ...
    'conditions does not match the\ndifferential order(s) appearing in the ' ...
    'problem.'])

% Need to sort the results of INITVALS above. This is because we can't guarantee
% that the LBC or RBC were imposed in the order of ascending variables, and
% ascending diffOrder. The SORTCONDITIONS() method of the TREEVAR class
% evaluates the conditions with TREEVAR inputs, which gives it enough
% information to be able to sort them in the correct order.
if ( isIVP )
    idx = treeVar.sortConditions(N.lbc, N.domain, diffOrders);
else
    idx = treeVar.sortConditions(N.rbc, N.domain, diffOrders);
end

% Sort the results from above:
initVals = initVals(idx);

% Create an ODESET struct for specifying tolerance:
opts = odeset('absTol', pref.ivpAbsTol, 'relTol', pref.ivpRelTol);

% What happiness check do we want to use for the IVP?
opts.happinessCheck = pref.happinessCheck;

% Do we want to restart the solver at breakpoints?
opts.restartSolver = pref.ivpRestartSolver;

% Break out of solver if solution exceeds maximum norm?
maxNorm = N.maxnorm;

% TODO: Should just need this in case opts.Events is not Inf. Move to end of
% file?
    function [position,isterminal,direction] = applyEventsFcn(t,y)
        position = abs(y(varIndex))-maxNorm(:); % The value that we want to be zero
        isterminal = 1+0*maxNorm;  % Halt integration in call cases
        direction = 0;   % The zero can be approached from either direction
    end

if ( ~isempty(maxNorm) && ~all(isinf(maxNorm) ))
    opts.Events = @applyEventsFcn;
end

% Solve!
[t, y]= solver(anonFun, odeDom, initVals, opts);

% To fit in with CHEBOP backslash, just return the functions, not their
% derivatives as well, as is normally done in the MATLAB ODE solvers. Above, we
% determined and stored in VARINDEX what columns of Y will correspond to the
% functions we want to return.
y = y(:, varIndex);

% Find how many dependent functions the solution consists of, i.e., the number
% of unknown variables in the case of coupled systems.
numColY = numColumns(y);

% If we were solving a coupled system, return the solution as a CHEBMATRIX to be
% consistent with standard CHEBOP backslash:
if (numColY > 1)
    y = chebmatrix(y);
    % Transpose that doesn't convert row chebfuns to column ones:
    y.blocks = reshape(y.blocks, numColY, 1); 
end

% Return what solver we were using
info.solver = solver;

% Do we want to display information about the solution process?
if ( ~strcmpi(pref.display, 'off') )
    displayInfo(y, isIVP)
end

% Return useful information about the solution:
info.solver = solver;

% Return a CHEBFUN rather than a CHEBMATRIX for scalar problems:
if ( ~isa(y, 'chebmatrix') )
    varargout{1} = y;
    varargout{2} = info;
elseif ( nargout == 1 )
    varargout{1} = y;
elseif ( nargout == size(y, 1) )
    [varargout{1:nargout}] = deal(y);
elseif ( nargout == size(y, 1) + 1 )
    [varargout{1:nargout - 1}] = deal(y);
    varargout{nargout} = info;
else
    error('CHEBFUN:CHEBOP:solveivp:numberOfOutputs', ...
        'Incorrect number of outputs.');
end

end
