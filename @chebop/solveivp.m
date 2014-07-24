function y = solveivp(N, varargin)

% Are we dealing with a system?
isSystem = ( nargin(N.op) <= 2 );

% Convert to first order format
if ( isSystem )
    anonFun = treeVar.toFirstOrder(N.op, N.domain);
    varIndex = 1;
else
    [anonFun, varIndex] = treeVar.toFirstOrderSystem(N.op, N.domain);
end
%% Obtain information about the initial conditions.
% Begin by evaluation N.lbc with a zero chebfun to pick up the desired values:

% Create a zero CHEBFUN for evaluating N.LBC or N.RBC. This gives us the correct
% values we need for passing as initial/final conditions to the ODE solver.
cheb0 = chebfun(@(x) 0*x, N.domain);

% Evaluate N.LBC or N.RBC:
if ( ~isempty(N.lbc) )
    cheb0 = repmat({cheb0}, nargin(N.lbc), 1);
    disp('Initial value problem detected.')
    bcEvalFun = N.lbc(cheb0{:});
    evalPoint = N.domain(1);
    odeDom = N.domain;
else
    cheb0 = repmat({cheb0}, nargin(N.rbc), 1);
    disp('Final value problem detected.')
    bcEvalFun = N.rbc(cheb0{:});
    evalPoint = N.domain(end);
    odeDom = fliplr(N.domain);
end

% Should sort the results. This requires evaluating N.LBC/RBC with a treeVar,
% looking at the IDs and the diffOrders.
if ( isa(bcEvalFun, 'chebfun') )
    initVals = -bcEvalFun(evalPoint);
else
    initVals = -cellfun(@feval, bcEvalFun.blocks, ...
        repmat({evalPoint}, length(bcEvalFun.blocks), 1) );
end

% Create an ODESET struct for specifying tolerance:
opts = odeset('absTol',1e-12,'relTol',1e-12);

[t, y]=chebfun.ode113(anonFun, odeDom, initVals, opts);
% To fit in with CHEBOP backslash, just return the functions, not their
% derivatives as well, as is normally done in the MATLAB ode solvers. Above, we
% determined and stored in VARINDEX what columns of Y will correspond to the
% functions.
y = y(:, varIndex);

% Return the solution as a CHEBMATRIX to be consistent with standard CHEBOP \:
if (numColumns(y) > 1)
    y = chebmatrix(y)';
end

end