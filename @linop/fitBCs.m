function u0 = fitBCs(L, prefs)
%FITBCS    Find a low-order polynomial which satisfies the BCs of a LINOP.
%   U0 = FITBCS(L) Returns a CHEBMATRIX which will satisfy the BCs
%        and other conditions of the linop L.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%  Determines how big discretization we want for each component. We want to
%  obtain the lowest degree polynomial which satisfies all boundary conditions.
%  In case of a system, this might mean that this won't be governed by
%  individual diffOrders, but rather, in how many BCs an unknown function
%  appears (e.g. u'+v = 0, u-v' = 0, u(-1) = u(1) = 0). This information can be
%  obtained from the iszero information of the linearised BCs in the linop L.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fix the discretization to @chebcolloc2, so that we can be sure we always get 
% the same initial guesses for Newton iteration.
prefs.discretization = @chebcolloc2;

% Store the total number of interior breakpoints
dom = L.domain;
numInts = length(dom) - 1;

% If we have subintervals, but no continuity conditions were passed, we need to
% create them.
if ( numInts > 1 && isempty(L.continuity) )
     L = deriveContinuity(L, dom);
end

% Number of variables appearing in the problem
numVar = size(L, 2);
isFun = isFunVariable(L);

% We will construct a low-degree polynomial for each unknown function in the
% problem. The degree of the polynomial depends on the number of conditions that
% the unknown function appears in, find that degree!
polyDegree = zeros(1, numVar);
if ( ~isempty(L.constraint.functional) )
    polyDegree = polyDegree + sum(~iszero(L.constraint.functional), 1);
end
% If continuity conditions are to be enforced, we need the polynomials to be at
% least the difforder of the variable it represents:
if ( ~isempty(L.continuity.functional) )
    polyDegree = max(polyDegree, max(L.diffOrder, [], 1)-1);
end
% Ensure that scalars only have length 1:
polyDegree(~isFun) = 1;
% And that each variable has degree at least 1:
polyDegree = max(polyDegree, 1);

% Create a discretization of the linear BCs:
discType = prefs.discretization;

% Set a size zero discretization. Actual size is controlled by dimAdjust.
dim = zeros(1, numInts);

% Initialize the discretization
B = 0;

% As a safeguard, we try finer discretizations until we have a discretized
% operator of a sufficient rank.
dimCounter = 0;
dimCounterMax = 5;
while ( rank(B) < size(B, 1) && dimCounter < dimCounterMax )
    
    % Create a discretization:
    disc = discType(L, dim);
    % Set its dimAdjust to required degree:
    disc.dimAdjust = polyDegree;

    % Obtain constraints and create the discrete (matrix) version of the BCs and
    % rhs values:
    B = getConstraints(disc);
    b = [];
    if ( ~isempty(L.constraint) )
        b = [ L.constraint.values ; b ];
    end
    if ( ~isempty(L.continuity) )
        b = [ L.continuity.values ; b ];
    end
    
    % Remove trivial rows:
    idx = ~any(B, 2);
    B(idx,:) = [];
    b(idx) = [];

    % Try increasing the discretization if we were not successful in getting a
    % full-rank B.
    dim = dim + 1;
    dimCounter = dimCounter + 1;

end


if ( dimCounter == dimCounterMax )
    % We failed. Return a zero initial guess.
    
    zeroFun = chebfun(0, dom);
    % Convert to a chebmatrix of correct dimensions
    u0 = cell(numVar, 1);
    for k = 1:numVar
        if ( isFun(k) )
            u0{k} = zeroFun;
        else
            u0{k} = 0;
        end
    end
    u0 = chebmatrix(u0);
   
    warning('CHEBFUN:LINOP:fitBCs:failure', ...
        'Unable to construct a suitable initial guess. Using a zero guess.')
    return
    
end

% Solve for the discrete values of the initial guess:
u0disc = B\(-b); % TODO: Why must b be negated?

% Chop u0disc into pieces.
u0disc = mypartition(disc, u0disc, polyDegree);

% Convert to a cell-array of CHEBFUN objects:
u0 = cell(numel(u0disc),1);
for k = 1:numel(u0)
    tmpDisc = disc;
    if ( isFun(k) )
        tmpDisc.dimension = tmpDisc.dimension + polyDegree(k);
    end
    u0{k} = toFunctionIn(tmpDisc, u0disc{k});
end

% Convert the cell-array of CHEBFUN objects to a CHEBMATRIX
u0 = chebmatrix(u0);

end

function u = mypartition(disc, values, dimAdjust)
%CHEBDISCRETIZATION.PARTITION   Partition values to appropriate variables.
%   U = CHEBDISCRETIZATION.PARTITION(DISC, VALUES) will, given a vector or
%   matrix (columnwise) VALUES of values corresponding to all the discretized
%   variables and scalars in a system DISC, convert to a cell-valued partition
%   of individual variables in the system. I.e., deduce the variable boundaries
%   within the discretization.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Document dimAdjust.

% Which variables are functions (as opposed to scalars)?
isFun = isFunVariable(disc.source);

if ( nargin < 3 )
    dimAdjust = zeros(size(isFun));
end

% Allocate the discretization size to each function variable.
m = ones(size(isFun));
m(isFun) = sum(disc.dimension) + dimAdjust(isFun)*disc.numIntervals;

% Do the partition.
u = mat2cell(values, m, size(values, 2));

end
