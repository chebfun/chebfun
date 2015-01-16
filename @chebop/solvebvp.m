function [u, info] = solvebvp(N, rhs, varargin)
%SOLVEBVP   Solve a linear or nonlinear CHEBOP BVP system.
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
%       ISLINEAR: A vector with four entries containing linearity information
%           for N. More specifically,
%               ISLINEAR(1) = 1 if N.OP is linear
%               ISLINEAR(2) = 1 if N.LBC is linear
%               ISLINEAR(3) = 1 if N.RBC is linear
%               ISLINEAR(4) = 1 if N.BC is linear
%           Otherwise, the corresponding element of ISLINEAR is equal to 0.
%
%   For linear problems, INFO further contains the field
%       ERROR:    The residual of the differential equation.
%
%   For nonlinear problems, INFO further contains the fields
%       NORMDELTA:  A vector of the norm of the Newton updates.
%       ERROR:      An error estimate for the convergence of the Newton
%                   iteration.
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
%   CHEBOP/SOLVEBVPNONLINEAR, CHEBOP/SOLVEIVP, LINOP/MLDIVIDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Developer note:
%   U = SOLVEBVP(N, RHS, PREF, DISPLAYINFO) allows passing in a function handle
%   to a displaying method that is called during the damped Newton iteration.
%   This allows separating the displaying process for regular CHEBOP use and
%   CHEBGUI. See chebop/displayInfo() and chebgui/displayInfo() for more
%   details.

% Parse inputs:
[pref, isPrefGiven, displayInfo] = parseInputs(N, varargin{:});

% Find out how many variables N operates on:
nVars = numVars(N);

% Support single input argument for autonomous scalar problems:
if ( nargin(N) == 1 )
    N.op = @(x, u) N.op(u);
end

% Store the domain we're working with.
dom = N.domain;

% Create an initial guess if none is passed.
if ( isempty(N.init) )
    % Initialise a zero CHEBFUN:
    zeroFun = chebfun(0, dom); 
    % Convert to a CHEBMATRIX of correct dimensions:
    u0 = cell(nVars, 1);
    for k = 1:nVars
        u0{k} = zeroFun;
    end
    u0 = chebmatrix(u0);
else
    % Get the initial guess.
    u0 = N.init; 
end

% Initialise the independent variable:
x = chebfun(@(x) x, dom);

% Linearize and attach preferences.
[L, residual, isLinear] = linearize(N, u0, x);

warnState = warning();
[ignored, lastwarnID] = lastwarn(); %#ok<ASGLU>
if ( strcmp(lastwarnID, 'CHEBFUN:CHEBOP:linearize:bcConcat') )
    warning('off', 'CHEBFUN:CHEBOP:linearize:bcConcat');
end

% Check the size of the residual (the output the dimensions of the CHEBOP).
[numRow, numCol] = size(residual);

% If the RHS passed is numerical, cast it to a CHEBMATRIX of the appropriate
% size before continuing:
if ( isnumeric(rhs) )
    % Check whether dimensions match:
    if ( ~(all(size(rhs) == [numRow, numCol])) &&  (max(size(rhs)) > 1) )
        if ( all(size(rhs) == [numCol, numRow]) )
            warning('CHEBFUN:CHEBOP:solvebvp:vertcat1', ...
                'Please concatenate the right-hand side of the BVP vertically. Transposing.')
            rhs = rhs.';
        else
            error('CHEBFUN:CHEBOP:solvebvp:rhs', ...
                'The right-hand side does not match the output dimensions of the operator.');
        end
    end
    
    % Convert the rhs to a CHEBMATRIX.
    rhs = rhs + 0*residual;
    
elseif ( isa(rhs, 'chebfun') && size(rhs, 2) > 1 )
    rhs = chebmatrix(mat2cell(rhs).');
    warning('CHEBFUN:CHEBOP:solvebvp:vertcat2', ...
        'Please use vertical concatenation for the right-side data.')
end

% Do the same for the initial guess:
if ( isnumeric(u0) )
    % Check whether dimensions match:
    if ( ~all(size(u0) == [numRow, numCol]) )
        if ( all(size(u0) == [numCol, numRow]) )
            warning('CHEBFUN:CHEBOP:solvebvp:vertcat3', ...
                ['Please concatenate the initial guess of the solution for '...
                'the BVP vertically. Transposing.']);
            u0 = u0.';
        else
            error('CHEBFUN:CHEBOP:solvebvp:init', ...
                'Initial guess does not match output dimensions of operator.');
        end
    end
    
    % Convert the initial guess to a CHEBMATRIX.
    u0 = u0 + 0*residual;
end

% Determine the discretization.
pref = determineDiscretization(N, L, isPrefGiven, pref);
disc = pref.discretization();

% Determine the TECH used by the discretization.
tech = disc.returnTech();
techUsed = tech();

% If the dicretization uses periodic functions, then clear the boundary
% conditions (if we're using periodic basis functions, the boundary conditions
% will be satisfied by construction). Also, ensure that u0 is of correct
% discretization, and convert it to a CHEBMATRIX if necessary.
if ( isPeriodicTech(techUsed) )
    % Clear the boundary conditions.
    [N, L] = clearPeriodicBCs(N, L);
    % Do the conversion.
    if ( isa(u0, 'chebfun') )
        u0 = chebmatrix(changeTech(u0, tech));
    elseif ( isa(u0, 'chebmatrix') )
        u0 = changeTech(u0, tech);
    end
end

% Solve:
if ( all(isLinear) )
    % Call solver method for linear problems.
    [u, info] = N.solvebvpLinear(L, rhs - residual, N.init, pref, displayInfo);
else
    % [TODO]: Switch between residual and error oriented Newton methods.
    
    % Create initial guess which satisfies the linearised boundary conditions:
    if ( isempty(N.init) )
        
        if ( ~isPeriodicTech(techUsed) )
            % Find a new initial guess that satisfies the BCs of L.
            % If we are using TRIGCOLLOC, we don't need to do that because 
            % the zero CHEBFUN is periodic.
            u0 = fitBCs(L, pref);
        end
        
        % Linearize about the new initial guess. If we are working with
        % parameter dependent problems, and did not get an initial condition
        % passed, we might have to cast some components in the CHEBMATRIX U0
        % from a CHEBFUN to a scalar. Hence, call LINEARIZE() with four outputs.
        [L, residual, isLinear, u0] = linearize(N, u0, x);
    end
    
    % If using a periodic TECH, ensure that rhs is of correct 
    % discretization, and convert it to a CHEBMATRIX if necessary.
    if ( isPeriodicTech(techUsed) )
        if ( isa(rhs, 'chebfun') )
            rhs = chebmatrix(changeTech(rhs, tech));
        elseif ( isa(rhs, 'chebmatrix') )
          rhs = changeTech(rhs, tech);
        end
    end

    % Call solver method for nonlinear problems.
    [u, info] = solvebvpNonlinear(N, rhs, L, u0, residual, pref, displayInfo);
    
end

% Revert warning state:
warning(warnState);

% Return a CHEBFUN rather than a CHEBMATRIX for scalar problems:
if ( all(size(u) == [1 1]) )
    u = u{1};
end

% Return the linearity information as well:
info.isLinear = isLinear;

end

function [pref, isPrefGiven, displayInfo] = parseInputs(N, varargin)
%PARSEINPUTS   Parse the input arguments to SOLVEBVP.

% Initialise the outputs:
pref = [];
displayInfo = [];

% Loop over varargin:
while ( ~isempty(varargin) )
    if ( ischar(varargin{1}) )
        warning('CHEBFUN:CHEBOP:solvebvp:stringInput', ...
            'String inputs to SOLVEBVP are deprecated.');
        varargin(1) = [];
    elseif ( isa(varargin{1}, 'cheboppref') )
        pref = varargin{1};
        varargin(1) = [];
        isPrefGiven = 1;
    elseif ( isa(varargin{1}, 'function_handle') )
        displayInfo = varargin{1};
        varargin(1) = [];
    else
        error('CHEBFUN:CHEBOP:solvebvp', ...
            'Unknown input of type %s to SOLVEBVP.', class(varargin{1}));
    end
end

% No preferences passed; use the current chebopprefs:
if ( isempty(pref) )
    pref = cheboppref();
    isPrefGiven = 0;
end

% If no DISPLAYINFO function handle passed, use the default CHEBOP one.
if ( isempty(displayInfo) )
    displayInfo = @chebop.displayInfo;
end

end
