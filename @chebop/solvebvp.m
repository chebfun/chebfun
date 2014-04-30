function [u, info] = solvebvp(N, rhs, pref, displayInfo)
%SOLVEBVP  Solve CHEBOP BVP system.
%
%   Note that CHEBOP requires the RHS of coupled systems to match the
%   system, even for scalars right-hand sides, e.g.,
%       N = chebop(@(x, u, v) [diff(u) + v ; u + diff(v)]);
%       N.bc = @(x, u, v) [u(-1) ; v(1)];
%       uv = solvebvp(N, 0);
%   is not an accepted syntax.
%
% See also: CHEBOP/MLDIVIDE.

% TODO: This is a user-facing function. It requires much better documentation.
            
% No preferences passed, use the current chebopprefs
if ( nargin < 3 )
    pref = cheboppref;
end

% If no DISPLAYINFO function handle passed, use the default CHEBOP one.
if ( nargin < 4 )
    displayInfo = @N.displayInfo;
end

% Support single input argument for autonomous scalar problems:
if ( nargin(N) == 1 )
    N.op = @(x, u) N.op(u);
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
    % Ensure that N.init is a CHEBMATRIX, not a CHEBFUN
    if ( isa(u0, 'chebfun') )
        u0 = chebmatrix(u0);
    end
end

% Initialise the independent variable:
x = chebfun(@(x) x, dom);

% Linearize
[L, residual, isLinear] = linearize(N, u0, x);

% Attach preferences to the LINOP:
L.prefs = pref;

% If the RHS passed is numerical, cast it to a CHEBMATRIX of appropriate size
% before continuing
if ( isnumeric(rhs) )
    rhs = N.convertToRHS(rhs, residual);
end

% Solve:
if ( all(isLinear) )
    % Call solver method for linear problems.
    [u, info] = N.solvebvpLinear(L, rhs, residual, displayInfo, pref);
    
else
    % TODO: Switch between residual and error oriented Newton methods
    
    % Create initial guess which satisfies the linearised boundary conditions:
    if ( isempty(N.init) )
        u0 = fitBCs(L);
        % Linearize about the new initial guess:
        [L, residual, isLinear] = linearize(N, u0, x);
    end

    % Call solver method for nonlinear problems.
    [u, info] = solvebvpNonlinear(N, rhs, L, u0, residual, pref, displayInfo);
    
end

% Return a CHEBFUN rather than a CHEBMATRIX for scalar problems:
if ( all(size(u) == [1 1]) )
    u = u{1};
end

% Return the linearity information as well
info.isLinear = isLinear;

end
