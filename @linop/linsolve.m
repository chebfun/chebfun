function [u, disc] = linsolve(L, f, discType)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( nargin < 3 )
    % TODO: Get from a (global?) preference?
    discType = L.discretizer;
end

isFun = isFunVariable(L); 

if isa( f, 'chebfun' )
    f = chebmatrix(f);
end

%% Set up the discretisation:
% Set the allowed discretisation lengths: (TODO: A preference?)
dimVals = floor(2.^[5 6 7 8 8.5 9 9.5 10 10.5 11]);

if ( isa(discType, 'function_handle') )
    % Create a discretization object
    disc = discType(L);  
    
    % Merge domains of the operator and the rhs:
    disc = mergeDomains(disc,f.domain); 
        
    % Update the discretization dimension on unhappy pieces:
    disc.dimension = repmat(dimVals(1), 1, numel(disc.domain)-1);
    dimVals(1) = [];
else
    % A discretisation is given. The idea is that it probably has an LU
    % factorization already attached, so try to use it first. Caller beware!
    disc = discType;
        
    % Initialise dimVals. Try the given size first, then iterate on the rest. 
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

