function u = solvebvp(N, rhs, pref)

% No preferences passed, use the current chebopprefs
if nargin < 3
    pref = cheboppref;
end

% NUMVARS indicate how many unknown function we seek.
numVars = nargin(N.op) - 1;

% Store the domain we're working with.
dom = N.domain;

% Initialise a zero ADCHEBFUN:
zeroFun = chebfun(0, dom);
u0 = cell(numVars, 1);
for k = 1:numVars
    u0{k} = zeroFun;
end
u0 = chebmatrix(u0);

% Initialise the dependent variable:
x = chebfun(@(x) x, dom);

% Linearise
[L, affine, isLinear] = linearise(N, x, u0);

% Solve:
if ( all(isLinear) )
    % TODO: Should this be calling the linop solve method directly (rather than
    % use backslash)?
    
    % Set the preferred discretization for the linop
    L.discretizationType = pref.discretizationType;
    
    % Solve the linear problem
    u = L\(rhs - affine);
    
    % TODO: Return residual as well?
else
    % Call solver method for nonlinear problems.
    % TODO: Swith between residual and error oriented Newton methods
    u = solvebvpNonlinear(N, rhs, pref);
end

% If we were solving a scalar problem, return a chebfun rather than a
% chebmatrix.
if ( all(size(u) == [1 1]) )
    u = u{1};
end

end