function [y, info] = solveivp(N, rhs, pref, varargin)

% Check inputs
if ( nargin < 2 )
    rhs = 0;
end

if ( nargin < 3 ) 
    pref = cheboppref;
end

if ( nargin < 4 )
    displayInfo = @chebop.displayIVPinfo;
else
    displayInfo = varargin{1};
end

% Check for unbounded domains:
if ( ~all(isfinite(N.domain)) )
    error('CHEBFUN:CHEBOP:solveivp:infDom', ...
        'Solving IVPs on unbounded intervals is not supported.');
end

% If pref.ivpSolver is set to a global method, we really should be calling
% CHEBOP/SOLVEBVP():
if ( isempty(strfind(func2str(pref.ivpSolver), 'chebfun.ode')) )
    y = solvebvp(N, rhs, pref, varargin{:});
    return
end

% Are we dealing with a system?
isSystem = ( nargin(N.op) <= 2 );

% Ensure RHS is a CHEBMATRIX
if ( ~isa(rhs, 'chebmatrix') )
    rhs = chebmatrix(rhs);
end

% Convert to first order format
if ( isSystem )
    [anonFun, varIndex, problemDom] = treeVar.toFirstOrder(N.op, rhs, N.domain);
else
    [anonFun, varIndex, problemDom] = treeVar.toFirstOrderSystem(N.op, rhs, N.domain);
end

% Join all breakpoints
dom = union(N.domain, problemDom);

%% Obtain information about the initial conditions.
% Begin by evaluating N.lbc with a zero chebfun to pick up the desired values:

% Create a zero CHEBFUN for evaluating N.LBC or N.RBC. This gives us the correct
% values we need for passing as initial/final conditions to the ODE solver.
cheb0 = chebfun(@(x) 0*x, dom);

% Evaluate N.LBC or N.RBC:
if ( ~isempty(N.lbc) )
    cheb0 = repmat({cheb0}, nargin(N.lbc), 1);
    bcEvalFun = N.lbc(cheb0{:});
    evalPoint = dom(1);
    odeDom = dom;
    % Store that we're solving an IVP.
    isIVP = 1;
else
    cheb0 = repmat({cheb0}, nargin(N.rbc), 1);
    bcEvalFun = N.rbc(cheb0{:});
    evalPoint = dom(end);
    % Flip the time domain, so that the problem will be solved as a final-value
    % problem.
    odeDom = fliplr(dom);
    % Store that we're solving an FVP.
    isIVP = 0;
end

% Evaluating N.LBC or N.RBC above gives CHEBFUN results. Evaluate the CHEBFUNs
% at the appropriate endpoint to obtain the scalar value for the initial/final
% condition.
if ( isa(bcEvalFun, 'chebfun') )
    initVals = -bcEvalFun(evalPoint);
else
    initVals = -cellfun(@feval, bcEvalFun.blocks, ...
        repmat({evalPoint}, length(bcEvalFun.blocks), 1) );
end

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

[t, y]= solver(anonFun, odeDom, initVals, opts);
% To fit in with CHEBOP backslash, just return the functions, not their
% derivatives as well, as is normally done in the MATLAB ode solvers. Above, we
% determined and stored in VARINDEX what columns of Y will correspond to the
% functions.
y = y(:, varIndex);

numColY = numColumns(y);
% Return the solution as a CHEBMATRIX to be consistent with standard CHEBOP \:
if (numColY > 1)
    y = chebmatrix(y);
    % Transpose that doesn't convert row chebfuns to column ones:
    y.blocks = reshape(y.blocks, numColY, 1); 
end

if ( ~strcmpi(pref.display, 'off') )
    % Do we want to display information about the solution process?
    displayInfo(y, isIVP)
end

end