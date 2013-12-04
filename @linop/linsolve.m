function [u, disc] = linsolve(L, f, discType)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( nargin < 3 )
    % TODO: Get from a (global?) preference?
    discType = L.discretizer;
end

isFun = isFunVariable(L); 

%% Set up the discretisation:
if ( isa(discType, 'function_handle') )
    % Create a discretization object
    disc = discType(L);  
    
    % MErge domains of the operator and the rhs:
    disc = mergeDomains(disc, L, f); % TODO: Why does this involve disc?
    
    % Set the allowed discretisation lengths: (TODO: A preference?)
    dimVals = floor(2.^[3 4 5 6 7 8 8.5 9 9.5 10 10.5 11]);
    
    % Update the discretistion dimension on unhappy pieces:
    disc.dimension = repmat(dimVals(1), 1, numel(disc.domain)-1);
    dimVals(1) = [];
else
    % A discretisation is given:
    disc = discType;
    
    % TODO: Check discretisation is valid for the given L and f!
    
    % Initialise dimVals;
    dimVals = max(disc.dimension);
end

if ( isempty(L.continuity) )
     % Apply continuity conditions:
     disc = deriveContinuity(disc);
end

% Initialise happiness:
numInt = disc.numIntervals;
isDone = false(1, numInt);

%% Loop over a finer and finer grid until happy:
for dim = dimVals

    % Discretize the operator (incl. constraints/continuity):
    A = disc.matrix();
    
    % Discretize the rhs (incl. constraints/continuity):
    b = disc.rhs(f);

    % Solve the linear system:
    [v, disc] = disc.mldivide(A, b);
    
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
%     u{k} = simplify(u{k}, epsLevel);
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
   
end

