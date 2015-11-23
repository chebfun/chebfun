function [u, disc] = linsolve(L, f, varargin)
%LINSOLVE  Solve a linear differential/integral equation.
%   Important: A CHEBOPPREF object PREFS has to be passed. When this method
%   is called via CHEBOP/MLDIVIDE, PREFS is inherited from the CHEBOP level.
%
%   U = LINSOLVE(L, F), or U = L\F, solves the linear system defined by L*U=F
%   for a LINOP L and CHEBMATRIX F. The result is a CHEBMATRIX.
%
%   An equivalent syntax to U = LINSOLVE(L, F) is U = L\F.
%
%   LINSOLVE(L,F,CDISC) uses the opDiscretization CDISC to solve the
%   problem. This can be used, for example, to introduce new breakpoints that
%   are not in the domain of either L or F.
%
%   LINSOLVE(...,PREFS) accepts a CHEBOPPREF to control the behavior of
%   the algorithms. If empty, defaults are used.
%
%   EXAMPLE:
%     d = [0,pi];
%     [Z,I,D] = linop.primitiveOperators(d);
%     A = linop( D^2 - I );
%     E = functionalBlock.eval(d);
%     A = addbc(A,E(0),0);
%     A = addbc(A,E(pi),1);
%     u = A \ chebfun('x',d);
%     plot(u{1})
%
% See also CHEBOPPREF, CHEBOP.MLDIVIDE.

%  Copyright 2015 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: The second output is a discretization that contains the
% stored LU factors that were used to find the final solution. If you pass
% that discretization back into another call, that factorization is tried
% first to save time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse input
prefs = [];    % no prefs given
disc = [];     % no discretization given
vscale = zeros(1, size(L, 2));
for j = 1:nargin-2
    item = varargin{j};
    if ( isa(item, 'cheboppref') )
        prefs = item;
    elseif ( isa(item,'opDiscretization') )
        disc = item;
    elseif ( isnumeric(item) )
        vscale = item(:)';
    else
        error('CHEBFUN:LINOP:linsolve:badInput', ...
            'Could not parse argument number %i.',j+2)
    end
end

% Check for unbounded domains:
if ( ~all(isfinite(L.domain)) )
    error('CHEBFUN:LINOP:linsolve:infDom', ...
        'Unbounded domains are not supported.');
end

% If RHS is a CHEBFUN or a DOUBLE, we need to convert it to CHEBMATRIX in
% order for the method to be able to work with it.
if ( isa( f, 'chebfun' ) )
    f = chebmatrix(f);
elseif ( isnumeric(f) )
    f = chebmatrix(mat2cell(f));
end

% Use a given discretization, or create one?
if ( isempty(disc) )
    % Construct the current globally set discretization:
    disc = prefs.discretization(L);
    % What values for the discretization do we want to consider?
    dimVals = disc.dimensionValues(prefs);
    % Update the domain if new breakpoints are needed:
    disc.domain = domain.merge(disc.domain, f.domain);
    % Update the dimensions to work with the correct number of breakpoints
    disc.dimension = repmat(dimVals(1), 1, numel(disc.domain) - 1);
    dimVals(1) = [];
else
    % We have to assume that the given L matches the discretization. Caller
    % beware!
    dimVals = disc.dimensionValues(prefs);
    dim1 = max(disc.dimension);
    dimVals = [ dim1, dimVals(dimVals > dim1) ];
end

% Derive automatic continuity conditions if none were given.
if ( isempty(L.continuity) )
    L = deriveContinuity(L, disc.domain);
    disc.source = L;
end

% Initialise happiness:
numInt = disc.numIntervals;
isDone = false(1, numInt);

%% Loop over a finer and finer grid until happy.
% We need to know which solution components to check for happiness:
isFun = isFunVariable(L);

for dim = [dimVals inf]

    % [TODO]: It's weird that the current value of dim is the _next_ disc size.
    % Discretize the operator (incl. constraints/continuity), unless there is a
    % currently valid factorization at hand.
    if ( isFactored(disc) )
        A = [];
        P = speye(disc.dimension*sum(isFun)+sum(~isFun));
    else
        [A, P] = matrix(disc);
        if ( size(A, 1) ~= size(A, 2) )
            warning('CHEBFUN:LINOP:linsolve:notSquare', ...
                ['Operator may not have the correct number of boundary ' ...
                'conditions.\n' ...
                'Matrix is not square. Problem may be ill-posed.']);
        end
    end
    
    % Discretize the RHS (incl. constraints/continuity):
    b = rhs(disc, f);
    
    % Solve the linear system:
    [v, disc] = mldivide(disc, A, b);
    
    % Project the solution:
    v = P*v;
    
    % [TODO]: We could test each variable at their input dimension, but then
    % each would be different and we would not be able to use the trick of
    % taking a linear combination. Instead we project and test convergence
    % at the size of the output dimension.
    
    % Convert the different components into cells:
    u = partition(disc, v);
    
    % Need a vector of vscales.
    if ( numel(vscale) == 1 ) 
        vscale = repmat(vscale, 1, length(isFun));
    end

    % Test the happiness of the function pieces:
    [isDone, cutoff, vscale] = ...
        testConvergence(disc, u(isFun), vscale(isFun), prefs);
    
    if ( all(isDone) || isinf(dim) )
        break
    else
        % Update the discretiztion dimension on unhappy pieces:
        disc.dimension(~isDone) = dim;
    end
    
end

if ( ~all(isDone) )
    warning('CHEBFUN:LINOP:linsolve:noConverge', ...
        'Linear system solution may not have converged.')
end

%% Tidy the solution for output:
% The variable u is a cell array with the different components of the solution.
% Because each function component may be piecewise defined, we will loop through
% one by one.
values = cat(2, u{isFun});
uOut = cell(size(values, 2), 1);
for k = 1:size(values, 2)
    v = disc.toFunctionOut(values(:,k),cutoff);
    uOut{k} = v;
end
u(isFun) = uOut;

% Convert to chebmatrix
u = chebmatrix(u);

end
