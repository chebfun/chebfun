function [u, disc] = linsolve(L, f, varargin)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

pref = L.prefs;
disc = [];
for j = 1:nargin-2
    item = varargin{j};
    if isa(item,'chebpref')
        %pref = chebpref(pref,item);
        error('Preferences must be set to the ''prefs'' property of the linop.')
    elseif isa(item,'chebDiscretization')
        disc = item;
    end
end

if isa( f, 'chebfun' )
    f = chebmatrix(f);
end

% Use a given discretization, or create one?
dimVals = pref.dimensionValues;
if isempty(disc)
    disc = pref.discretization(L);
    disc = mergeDomains(disc,f.domain); 
    disc.dimension = repmat(dimVals(1), 1, numel(disc.domain)-1);
    dimVals(1) = [];
else
    dim1 = max(disc.dimension);
    dimVals = [ dim1, dimVals(dimVals > dim1) ];
end

% Derive automatic continuity conditions if none were given.
if ( isempty(L.continuity) )
     disc.source = deriveContinuity(L);
end

% Initialise happiness:
numInt = disc.numIntervals;
isDone = false(1, numInt);

%% Loop over a finer and finer grid until happy:
isFun = isFunVariable(L); 
for dim = dimVals

    % Discretize the operator (incl. constraints/continuity), unless there is a
    % currently valid factorization at hand. 
    if ( ~isFactored(disc) )
        A = matrix(disc);
    else
        A = [];
    end
    
    % Discretize the rhs (incl. constraints/continuity):
    b = rhs(disc,f);
 
    % Solve the linear system:
    [v, disc] = mldivide(disc,A, b);
    
    % Convert the different components into cells
    u = partition(disc,v);
   
    % Test the happieness of the function pieces:
    [isDone, epsLevel] = testConvergence(disc, u(isFun));
    
    if ( all(isDone) )
        break
    else
        % Update the discretistion dimension on unhappy pieces:
        disc.dimension(~isDone) = dim;
    end
    
end

if ( ~all(isDone) )
    warning('LINOP:linsolve:NoConverge', ...
        'Linear system solution may not have converged.')
end

%% Tidy the solution for output:
% The variable u is a cell array with the different components of the solution.
% Because each function component may be piecewise defined, we will loop through
% one by one.
for k = find( isFun )
    u{k} = disc.toFunction(u{k}); 
    u{k} = simplify(u{k}, epsLevel);
end

% Convert to chebmatrix
u = chebmatrix(u);

end

%% Auxillary functions:

function u = myNum2Cell(data, dim, isFunVariable)
% Break discrete solution into chunks representing functions and scalars:

    m = ones(size(isFunVariable));
    m(isFunVariable) = sum(dim);
    u = mat2cell(data, m, 1);
    uFun = u(isFunVariable);
   
end
