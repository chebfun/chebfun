function [u, info] = solvebvp(N, rhs, pref, displayInfo)
%SOLVEBVP  Solve a linear or nonlinear CHEBOP BVP system.
%
%   U = SOLVEBVP(N, RHS), where N is a CHEBOP and RHS is a CHEBMATRIX, CHEBFUN
%   or a vector of doubles attempts to solve the BVP
%
%       N(U) = RHS + boundary conditions specified by N
%
%   Observe that U = SOLVEBVP(N, RHS) has the same effect as U = N\RHS, but this
%   method allows greater flexibility than CHEBOP backslash, as described below.
%
%   If successful, the solution returned, U, is a CHEBFUN if N specifies a
%   scalar problem, and a CHEBMATRIX if N specifies a coupled systems of
%   ordinary differential equations. If N specifies a linear operator, the BVP
%   is solved using a spectral or a pseudospectral method. If N specifies a
%   nonlinear operator, damped Newton iteration in function space is performed,
%   where each linear problem arising is solved via a spectral/pseudospectral
%   method.
%
%   U = SOLVEBVP(N, RHS, PREF) is the same as above, using the preferences
%   specified by the CHEBOPPREF variable PREF.
%
%   [U, INFO] = SOLVEBVP(N, RHS, PREF) is the same as above, but also returns
%   the MATLAB struct INFO, which contains useful information about the solution
%   process. The fields of INFO are as follows:
%       ERROR:    The residual of the differential equation.
%       ISLINEAR: A vector with four entries containing linearity information
%           for N. More specifically, 
%               ISLINEAR(1) = 1 if N.OP is linear
%               ISLINEAR(2) = 1 if N.LBC is linear
%               ISLINEAR(3) = 1 if N.RBC is linear
%               ISLINEAR(4) = 1 if N.BC is linear
%           Otherwise, the corresponding element of ISLINEAR is equal to 0.
%
%   TODO: INFO will have more fields once we move into nonlinear problems,
%   update the list accordingly.
%
%   Note that CHEBOP allows the RHS of coupled system of ODEs to be a scalar,
%   e.g., one can both call
%       N = chebop(@(x, u, v) [diff(u) + v ; u + diff(v)]);
%       N.bc = @(x, u, v) [u(-1) ; v(1)];
%       uv = solvebvp(N, 0);
%   and
%       uv = solvebvp(N, [0; 0]);
%
% See also: CHEBOP, CHEBOP/MLDIVIDE, CHEBOPPREF, CHEBOP/SOLVEBVPLINEAR, 
%   CHEBOP/SOLVEBVPNONLINEAR, LINOP/MLDIVIDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Developers note:
%   U = SOLVEBVP(N, RHS, PREF, DISPLAYINFO) allows passing in a function handle
%   to a displaying method that is called during the damped Newton iteration.
%   This allows separating the displaying process for regular CHEBOP use and
%   CHEBGUI. See chebop/displayInfo() and chebgui/displayInfo() for more
%   details.

% No preferences passed; use the current chebopprefs:
if ( nargin < 3 )
    pref = cheboppref();
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
    % Convert to a chebmatrix of correct dimensions:
    u0 = cell(numVars, 1);
    for k = 1:numVars
        u0{k} = zeroFun;
    end
    u0 = chebmatrix(u0);
else
    u0 = N.init;
    % Ensure that N.init is a CHEBMATRIX, not a CHEBFUN:
    if ( isa(u0, 'chebfun') )
        u0 = chebmatrix(u0);
    end
end

% Initialise the independent variable:
x = chebfun(@(x) x, dom);

% Linearize and attach preferences:
[L, residual, isLinear] = linearize(N, u0, x);

% Check the size of the residual (the output the dimensions of the CHEBOP).
[numRow, numCol] = size(residual);

% If the RHS passed is numerical, cast it to a CHEBMATRIX of the appropriate
% size before continuing:
if ( isnumeric(rhs) )
    % Check whether dimensions match:
    if ( ~(all(size(rhs) == [numRow, numCol])) &&  (max(size(rhs)) > 1) )
        if ( all(size(rhs) == [numCol, numRow]) )
            warning('CHEBFUN:CHEBOP:solvebvp', ...
                'Please concatenate RHS of the BVP vertically. Transposing.')
            rhs = rhs.';
        else
            error('CHEBFUN:CHEBOP:solvebvp:rhs', ...
               'RHS does not match output dimensions of operator.');
        end
    end
    
    % If we get here, we have something compatable, this is a simple way to
    % convert RHS to a CHEBMATRIX:
    rhs = rhs + 0*residual;
    
elseif ( isa(rhs, 'chebfun') && size(rhs, 2) > 1 )
    rhs = chebmatrix(mat2cell(rhs).');
    warning('CHEBFUN:CHEBOP:solvebvp:vertcat', ...
        'Please use vertical concatenation for RHS ')
end

% Do the same for the initial guess:
if ( isnumeric(u0) )
    % Check whether dimensions match:
    if ( ~all(size(u0) == [numRow, numCol]) )
        if ( all(size(u0) == [numCol, numRow]) )
            warning('CHEBFUN:CHEBOP:solvebvp', ...
                ['Please concatenate the initial guess of the solution for '...
                'the BVP vertically. Transposing.']);
            u0 = u0.';
        else
            error('CHEBFUN:CHEBOP:solvebvp:init', ...
                'Initial guess does not match output dimensions of operator.');
        end
    end
    
    % Convert the initial guess to a CHEBMATRIX
    u0 = u0 + 0*residual;
end

% Solve:
if ( all(isLinear) )
    % Call solver method for linear problems.
    [u, info] = N.solvebvpLinear(L, rhs - residual, pref, displayInfo);
    
else
    % TODO: Switch between residual and error oriented Newton methods.
    
    % Create initial guess which satisfies the linearised boundary conditions:
    if ( isempty(N.init) )
        u0 = fitBCs(L, pref);
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

% Return the linearity information as well:
info.isLinear = isLinear;

end