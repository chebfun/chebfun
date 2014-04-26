function [u, info] = solvebvp(N, rhs, pref, displayInfo)

% No preferences passed, use the current chebopprefs
if ( nargin < 3 )
    pref = cheboppref;
end

% If no DISPLAYINFO function handle passed, use the default CHEBOP one.
if ( nargin < 4 )
    displayInfo = @N.displayInfo;
end

% NUMVARS indicate how many unknown function we seek.
numVars = max(nargin(N) - 1, 1);

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
        u0 = chebmatrix(u0);
    end
end
% Initialise the dependent variable:
x = chebfun(@(x) x, dom);

% Linearize
[L, residual, isLinear] = linearize(N, x, u0);

% If the RHS passed is numerical, cast it to a CHEBMATRIX of appropriate size
% before continuing
if ( isnumeric(rhs) )
    rhs = N.convertToRHS(rhs, residual);
end

% Attach the preferences to the linop before continuing
L.prefs = pref;

% Solve:
if ( all(isLinear) )
    % Call solver method for linear problems.
    [u, info] = solvebvpLinear(N, L, rhs, residual, x, displayInfo, pref);
else
    % Call solver method for nonlinear problems.
    % TODO: Swith between residual and error oriented Newton methods
    
    % Create initial guess which satisfies the linearised boundary conditions:
    if ( isempty(N.init) )
        u0 = fitBCs(L);
        % Linearize about the new initial guess:
        [L, residual, isLinear] = linearize(N, x, u0);
    end

    [u, info] = solvebvpNonlinear(N, rhs, L, u0, residual, pref, displayInfo);
end

% If we were solving a scalar problem, return a CHEBFUN rather than a
% CHEBMATRIX.
if ( all(size(u) == [1 1]) )
    u = u{1};
end

% Return the linearity information as well
info.isLinear = isLinear;

end
