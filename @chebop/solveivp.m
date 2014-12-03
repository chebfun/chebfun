function [y, info] = solveivp(N, rhs, pref, varargin)
%SOLVEIVP    Solve an IVP by reforming it to a first order system.
%
%   U = SOLVEIVP(N, RHS), where N is a CHEBOP and RHS is a CHEBMATRIX, CHEBFUN
%   or a vector of doubles attempts to solve the IVP
%
%       N(U) = RHS + boundary conditions specified by N
%
%   Observe that U = SOLVEIVP(N, RHS) has the same effect as U = N\RHS, but this
%   method allows greater flexibility than CHEBOP backslash, as described below.
%
%   If successful, the solution returned, U, is a CHEBFUN if N specifies a
%   scalar problem, and a CHEBMATRIX if N specifies a coupled systems of
%   ordinary differential equations. This method solves both linear and
%   nonlinear problems be automatically converting them to a coupled first-order
%   system, which can then be solved using MATLAB's built in solvers.
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
%
%   Note that CHEBOP allows the RHS of coupled system of ODEs to be a scalar,
%   e.g., one can both call
%       N = chebop(@(x, u, v) [diff(u) + v ; u + diff(v)], [0 10]);
%       N.bc = @(x, u, v) [u(0) ; v(0)];
%       uv = solvebvp(N, 0);
%   and
%       uv = solvebvp(N, [0; 0]);
%
% See also: CHEBOP, CHEBOP/MLDIVIDE, CHEBOPPREF, CHEBOP/SOLVEBVP,
% CHEBFUN/ODE113, CHEBFUN/ODE15S, CHEBFUN/ODE45, CHEBFUN/ODESOL, TREEVAR. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
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

% If pref.ivpSolver is set to a global method, we really should be calling
% CHEBOP/SOLVEBVP():
if ( isempty(strfind(func2str(pref.ivpSolver), 'chebfun.ode')) )
    [y, info] = solvebvp(N, rhs, pref, varargin{:});
    info.solver = 'Global method';
    return
end

%% Convert to a first-order system

% We call the conversion methods of the TREEVAR class, the call depends on
% whether we're dealing with a system or not. We don't support integral
% equations, so listen out for warnings that the TREEVAR class throws if it
% encounters unsupported methods.
try
    [anonFun, varIndex, problemDom, coeffs, diffOrders] = ...
        treeVar.toFirstOrder(N.op, rhs, N.domain);
catch ME
    % Did we encounter an unsupported method? If so, try to solve it globally:
    if ( ~isempty(regexp(ME.identifier, 'CHEBFUN:TREEVAR:.+:notSupported', ...
            'once')) )
        [y, info] = solvebvp(N, rhs, pref, varargin{:});
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
    % systems:
    cheb0 = repmat({cheb0}, nargin(N.lbc), 1);
    % Evaluate the initial condition:
    bcEvalFun = N.lbc(cheb0{:});
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
    cheb0 = repmat({cheb0}, nargin(N.rbc), 1);
    % Evaluate the final condition:
    bcEvalFun = N.rbc(cheb0{:});
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
    idx = treeVar.sortConditions(N.lbc, N.domain);
else
    idx = treeVar.sortConditions(N.rbc, N.domain);
end

% Sort the results from above:
initVals = initVals(idx);

% Create an ODESET struct for specifying tolerance:
opts = odeset('absTol', pref.ivpAbsTol, 'relTol', pref.ivpRelTol);

% What solver do we want to use for the IVP?
solver = pref.ivpSolver;

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

end