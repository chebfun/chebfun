function [u, info] = solvebvp(N, rhs, pref)

% No preferences passed, use the current chebopprefs
if nargin < 3
    pref = cheboppref;
end

% NUMVARS indicate how many unknown function we seek.
numVars = nargin(N) - 1;

% Store the domain we're working with.
dom = N.domain;

% Create an initial guess if none is passed
if ( isempty(N.init) )
    % Initialise a zero CHEBFUN:
    zeroFun = chebfun(0, dom);
    % Convert to a chebmatrix of correct dimensions
    u0 = cell(numVars, 1);
    for k = 1:numVars
        u0{k} = zeroFun;
    end
    u0 = chebmatrix(u0);
else
    u0 = N.init;
    % TODO: This should be done automatically when N.init is assigned
    if ( isa(u0, 'chebfun') )
        u0 = chebmatrix({u0});
    end
end
% Initialise the dependent variable:
x = chebfun(@(x) x, dom);

% Linearise
[L, residual, isLinear] = linearise(N, x, u0);

% Solve:
if ( all(isLinear) )
    % TODO: Should this be calling the linop solve method directly (rather than
    % use backslash)?
    
    % Set the preferred discretization for the linop
    L.discretizationType = pref.discretizationType;
    
    % Solve the linear problem
    u = L\(rhs - residual);
    
    % TODO: Return residual as well?
else
    % Call solver method for nonlinear problems.
    % TODO: Swith between residual and error oriented Newton methods
    [u, info] = solvebvpNonlinear(N, rhs, L, u0, residual, pref);
end

% If we were solving a scalar problem, return a chebfun rather than a
% chebmatrix.
if ( all(size(u) == [1 1]) )
    u = u{1};
end

end