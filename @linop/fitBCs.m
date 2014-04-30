function u0 = fitBCs(L)
%FITBCS    Find a low-order polynomial which satisfies the BCs of a LINOP.
%   U0 = FITBCS(L) Returns a CHEBMATRIX which will satisfy the BCs
%        and other conditions of the linop L.

% Store the total number of interior breakpoints
dom = L.domain;
numInts = length(dom) - 1;

% Determines how big discretization we want for each component. We want to
% obtain the lowest degree polynomial which satisfies all boundary
% conditions. In case of a system, this might mean that this won't be
% governed by individual diffOrders, but rather, in how many BCs an unknown
% function appears (e.g. u'+v = 0, u-v' = 0, u(-1) = u(1) = 0). This
% information can be obtained from the iszero information of the linearised
% BCs in the linop L.

% If we have subintervals, but no continuity conditions were passed, we need to
% create them.
if ( numInts > 1 && isempty(L.continuity) )
     L = deriveContinuity(L, dom);
end

% Number of variables appearing in the problem
numVar = size(L, 2);

% We will construct a low-degree polynomial for each unknown function in the
% problem. The degree of the polynomial depends on the number of conditions that
% the unknown function appears in, find that degree!
polyDegree = zeros(1, numVar);
if ( ~isempty(L.constraint.functional) )
    polyDegree = polyDegree + sum(~iszero(L.constraint.functional), 1);
end

if ( ~isempty(L.continuity.functional) )
    polyDegree = polyDegree + sum(~iszero(L.continuity.functional), 1);
end
polyDegree = 2;
polyDegree = max(max(polyDegree), 1);
dim = repmat(polyDegree(1), 1, numInts);

% Create a discretization of the linear BCs:
discType = L.prefs.discretization;

B = 0;
j = 0;
while ( rank(B) < size(B, 1) && j < 5 )
    
    disc = discType(L, dim);
    disc.dimAdjust = zeros(size(disc.dimAdjust));

    % Create the discrete (matrix) version of the BCs and rhs values:
    B = getConstraints(disc);
    b = [];
    if ( ~isempty(L.constraint) )
        b = [ L.constraint.values ; b ];
    end
    if ( ~isempty(L.continuity) )
        b = [ L.continuity.values ; b ];
    end

    dim = dim + 1;
    j = j + 1;
end

if ( j == 5 )
    % We failed. Returns a zero initial guess.
    
    zeroFun = chebfun(0, dom);
    % Convert to a chebmatrix of correct dimensions
    u0 = cell(numVar, 1);
    for k = 1:numVar
        u0{k} = zeroFun;
    end
    u0 = chebmatrix(u0);
   
    warning('Unable to construct a suitable initial guess. Using a zero guess.')
    return
end

% Solve for the discrete values of the initial guess:
u0disc = B\(-b); % TODO: Why must b be negated?

u0disc = partition(disc, u0disc);

% Convert to a cell-array of CHEBFUN objects:
u0 = cell(numel(u0disc),1);
for k = 1:numel(u0)
    u0{k} = disc.toFunction(u0disc{k}, 2);
end

% Convert the cell-array of CHEBFUN objects to a CHEBMATRIX
u0 = chebmatrix(u0);

end