function [u, disc] = linsolve(L, f, varargin)
%LINSOLVE  Solve a linear differential/integral equation.
%   Important: While you can construct a LINOP and apply this method, the 
%   recommended procedure is to use CHEBOP.MLDIVIDE instead. 
%   
%   U = LINSOLVE(L, F), or U = L\F, solves the linear system defined by L*U=F
%   for a linop L and chebmatrix F. The result is a chebmatrix.
%
%   Parameters controlling the method of solution are found and set using
%   L.prefs. 
%
%   U = LINSOLVE(L, F, CDISC) uses the chebDiscretization CDISC to solve the
%   problem. This can be used, for example, to introduce new breakpoints that
%   are not in the domain of either L or F.
%
%   See also CHEBOPPREF, CHEBOP.MLDIVIDE.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note: The second output is a discretization that contains the
% stored LU factors that were used to find the final solution. If you pass
% that discretization back into another call, that factorization is tried
% first to save time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Store preferences of the input LINOP
pref = L.prefs;
disc = [];

% Parse input
for j = 1:nargin-2
    item = varargin{j};
    if isa(item,'chebpref')
        error('Preferences must be set to the ''prefs'' property of the linop.')
    elseif isa(item,'chebDiscretization')
        disc = item;
    end
end

% If RHS is a CHEBFUN, need to convert it to CHEBMATRIX in order for the method
% to be able to work with it.
if isa( f, 'chebfun' )
    f = chebmatrix(f);
end

% Use a given discretization, or create one?
dimVals = pref.dimensionValues;
if isempty(disc)
    disc = pref.discretization(L);
    % Update the domain if new breakpoints are needed
    disc.domain = chebfun.mergeDomains(disc.domain, f.domain); 
    % Update the dimensions to work with the correct number of breakpoints
    disc.dimension = repmat(dimVals(1), 1, numel(disc.domain) - 1);
    dimVals(1) = [];
else
    % We have to assume that the given L matches the discretization. Caller
    % beware!
    dim1 = max(disc.dimension);
    dimVals = [ dim1, dimVals(dimVals > dim1) ];
end

% Derive automatic continuity conditions if none were given.
if ( isempty(L.continuity) )
     L = deriveContinuity(L,disc.domain);
     disc.source = L;
end

% Initialise happiness:
numInt = disc.numIntervals;
isDone = false(1, numInt);

%% Loop over a finer and finer grid until happy:
% TODO: What implications does the isFun variable have? 
isFun = isFunVariable(L); 
for dim = dimVals

    % Discretize the operator (incl. constraints/continuity), unless there is a
    % currently valid factorization at hand. 
    if ( ~isFactored(disc) )
        A = matrix(disc);
    else
        A = [];
    end
    
    % Discretize the RHS (incl. constraints/continuity):
    b = rhs(disc, f);
 
    % Solve the linear system:
    [v, disc] = mldivide(disc, A, b);
    
    % Convert the different components into cells
    u = partition(disc, v);
   
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

