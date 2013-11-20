function [u, disc] = linsolve(L, f, discType)

if ( nargin < 3 )
    % TODO: Get from a (global?) preference?
    discType = L.discretization;
end

[rowSize, colSize] = blockSizes(L.operator);
isFunVariable = isinf(colSize(1, :));

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

    % Solve the liner system:
    [v, disc] = disc.mldivide(A, b);
    
    % Convert the array of data to a cell
    u = myNum2Cell(v, disc.dimension, isFunVariable);
    
    % Take a inear combination of the function variables:
    uLin = linearCombination(u(isFunVariable), disc.dimension);
   
    % Test the happieness of the function pieces:
    [isDone, epsLevel] = myTestConvergence(disc, uLin);
    
    if ( all(isDone) )
        break
    else
        % Update the discretistion dimension on unhappy pieces:
        disc.dimension(~isDone) = dim;
    end
    
end

if ( ~isDone )
    warning('LINOP:linsolve:NoConverge', ...
        'Linear system solution may not have converged.')
end

%% Tidy the solution for output:
% The variable u is a cell array with the different components of the solution.
% Because each function component may be piecewise defined, we will loop through
% one by one.
for k = find( isFunVariable )
    u{k} = disc.toFunction(u{k}); 
%     u{k} = simplify(u{k}, epsLevel);
end

if ( numel(u) > 1)
    % Conver the cell array to a CHEBMATRIX:
    u = chebmatrix(u);
else
    % Output a CHEBFUN for scalar equations:
    u = u{1};
end

end

%% Auxillary functions:

function u = myNum2Cell(data, dim, isFunVariable)
% Break discrete solution into chunks representing functions and scalars:

    m = ones(size(isFunVariable));
    m(isFunVariable) = sum(dim);
    u = mat2cell(data, m, 1);
    uFun = u(isFunVariable);
   
end

function uLin = linearCombination(uFun, dim)
% ULIN returns a cell containing an arbitrary linear combination of the function
% variables, which each entry coresponding to different intervals.
    if ( nargin < 2 )
        dim = length(uFun{1});
    end
    % Take an arbitrary linear combination:
    s = 1 ./ (3*(1:numel(uFun))).';
    linVals = cell2mat(uFun.')*s;
    uLin = mat2cell(linVals, dim, 1);    
end

function [isDone, epsLevel] = myTestConvergence(disc, v)
% Test convergence of the solution between each of the breakpoints.:

% Test happiness:
numInt = numel(disc.domain)-1;
isDone = true(1, numInt);
epsLevel = 0;
for i = 1:numInt
    [isDone(i), t2] = disc.testConvergence(v{i});
    epsLevel = max(epsLevel, t2);
end
    
end