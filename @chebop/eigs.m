function varargout = eigs(N, varargin)
%EIGS   Find selected eigenvalues and eigenfunctions of a linear CHEBOP.
%   D = EIGS(A) returns a vector of 6 eigenvalues of the linear CHEBOP A. EIGS
%   will attempt to return the eigenvalues corresponding to the least
%   oscillatory eigenfunctions. (This is unlike the built-in EIGS, which returns
%   the largest eigenvalues by default.) If A is not linear, an error is
%   returned.
%
%   [V, D] = EIGS(A) returns a diagonal 6x6 matrix D of A's least oscillatory
%   eigenvalues, and a quasimatrix V of the corresponding eigenfunctions.
%
%   EIGS(A, B) solves the generalized eigenproblem A*V = B*V*D, where B is
%   another chebop on the same domain.
%
%   EIGS(A, K) and EIGS(A, B, K) for an integer K > 0 find the K smoothest
%   eigenvalues.
%
%   EIGS(A, K, SIGMA) and EIGS(A, B, K, SIGMA) find K eigenvalues. If SIGMA is a
%   scalar, the eigenvalues found are the ones closest to SIGMA. Other
%   possibilities are 'LR' and 'SR' for the eigenvalues of largest and smallest
%   real part, and 'LM' (or Inf) and 'SM' for largest and smallest magnitude.
%   SIGMA must be chosen appropriately for the given operator; for example, 'LM'
%   for an unbounded operator will fail to converge!
%
%   EIGS(..., PREFS) accepts a CHEBOPPREF to control the behavior of the
%   algorithm. If empty, defaults are used.
%
%   Despite the syntax, this version of EIGS does not use iterative methods
%   as in the built-in EIGS for sparse matrices. Instead, it uses the
%   built-in EIG on dense matrices of increasing size, stopping when the 
%   targeted eigenfunctions appear to have converged, as determined by the
%   chebfun constructor.
%
% Example:
%   N = chebop(@(u) diff(u, 2), [0 pi], 'dirichlet');
%   [V, D] = eigs(N, 10);
%   format long, sqrt(-diag(D))  % integers, to 14 digits
%   plot(V) % scaled sine waves
%
% See also LINOP/EIGS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Did we get preferences passed?
if ( (nargin > 1) && isa(varargin{end}, 'cheboppref') )
    prefs = varargin{end};
    isPrefGiven = 1;
else
    prefs = cheboppref();
    isPrefGiven = 0;
end

% Tell CHEBOP/LINEARIZE() to stop if it detects nonlinearity:
linCheck = true; 

% Linearize, thereby obtaining linearity information, a LINOP, and an input of
% the correct dimensions to pass to N:
[L, ~, isLinear, u0] = linearize(N, N.init, [], linCheck);

% We need the entire operator (including BCs) to be linear:
assert(all(isLinear), 'CHEBFUN:CHEBOP:eigs:nonlinear', ...
    ['The input operator appears to be nonlinear.\n', ...
    'EIGS() supports only linear CHEBOP instances.']);

% Support for generalised problems:
if ( nargin > 1 && isa(varargin{1}, 'chebop') )
    % Tell CHEBOP/LINEARIZE() that we don't want it to try to reshape inputs
    % that it believes are parameters to doubles, rather than CHEBFUNs.
    paramReshape = false;
    
    % Linearise the second CHEBOP:
    [varargin{1}, ~, isLinear] = ...
        linearize(varargin{1}, u0, [], linCheck, paramReshape);

    % We need the entire operator (including BCs) to be linear:
    assert(all(isLinear), 'CHEBFUN:CHEBOP:eigs:nonlinear', ...
        ['The second input operator appears to be nonlinear.\n', ...
        'EIGS() supports only linear CHEBOP instances.']);
    
end

% Determine the discretization.
prefs = determineDiscretization(N, length(L.domain), prefs);

% Clear boundary conditions if the dicretization uses periodic functions (since
% if we're using periodic basis functions, the boundary conditions will be
% satisfied by construction).
disc = prefs.discretization();
tech = disc.returnTech();
if ( isPeriodicTech(tech()) )
    [~, L] = clearPeriodicBCs(N, L);
end

% Add the preferences in vargarin to pass them to LINOP/EIGS.
if ( isPrefGiven )
    % If a CHEBOPPREF was passed to the method, it will have been at the last
    % position of varargin, indexed at nargin-1. Overwrite it with the current
    % PREFS, as the discretization might have changed in the periodic case:
    varargin{nargin-1} = prefs;
else
    % Otherwise, add the PREFS to VARARGIN, so that it can be passed to the call
    % to LINOP/EIGS below.
    varargin{nargin} = prefs;
end


% Call LINOP/EIGS.
[varargout{1:nargout}] = eigs(L, varargin{:});

% Return a CHEBFUN rather than a CHEBMATRIX for scalar problems:
if ( nargout > 1 && isa(varargout{1}, 'chebmatrix') )
    u = varargout{1};
    if ( all(size(u) == [1 1]) )
        u = u{1};
    end
    varargout{1} = u;
end

end
