function varargout = svds(N, varargin)
%SVDS   Find selected singularvalues and singularfunctions of a linear CHEBOP.
%   S = SVDS(A) returns a vector of 6 singularvalues of the linear CHEBOP A.
%   SVDS will attempt to return the singularvalues corresponding to the least
%   oscillatory singularfunctions. (This is unlike the built-in SVDS, which
%   returns the largest singularvalues by default.) If A is not linear, an
%   error is returned.
%
%   [U, S, V] = SVDS(A) returns a diagonal 6x6 matrix S of A's least
%   oscillatory singularvalues, and a quasimatrices U and V  of the
%   corresponding left and right singularfunctions.
%
%   [...] = SVDS(A, K) for an integer K > 0 find the K smoothest
%   singularvalues.
%
%   SVDS(..., PREFS) accepts a CHEBOPPREF to control the behavior of the
%   algorithm. If empty, defaults are used.
%
%   Despite the syntax, this version of SVDS does not use iterative methods
%   as in the built-in SVDS for sparse matrices. Instead, it uses the
%   built-in SVD on dense matrices of increasing size, stopping when the 
%   targeted singularfunctions appear to have converged, as determined by the
%   chebfun constructor.
%
% Example:
%   N = chebop(@(u) diff(u), [0 pi]);
%   S = svds(N, 10), format long  % integers, to 14 digits
%
% See also LINOP/SVDS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Did we get preferences passed?
if ( (nargin > 1) && isa(varargin{end}, 'cheboppref') )
    prefs = varargin{end};
    isPrefGiven = 1;
else
    prefs = cheboppref();
    isPrefGiven = 0;
end

% check the boundary conditions
bcType = getBCType(N);
if ( ~strcmp(bcType,'bvp') && ~strcmp(bcType,'periodic') )
    error('CHEBFUN:CHEBOP:svds:bcs', ...
    'SVDS only supports periodic or endpoint boundary conditions.');
end

% Tell CHEBOP/LINEARIZE() to stop if it detects nonlinearity:
linCheck = true; 

% Linearize, thereby obtaining linearity information, a LINOP, and an input of
% the correct dimensions to pass to N:
[L, ~, isLinear, u0] = linearize(N, N.init, [], linCheck);

% We need the entire operator (including BCs) to be linear:
assert(all(isLinear), 'CHEBFUN:CHEBOP:svds:nonlinear', ...
    ['The input operator appears to be nonlinear.\n', ...
    'SVDS() supports only linear CHEBOP instances.']);

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

% Add the bcType and preferences in vargarin to pass them to LINOP/SVDS.
if ( isPrefGiven )
    % If a CHEBOPPREF was passed to the method, it will have been at the last
    % position of varargin, indexed at nargin-1. Overwrite it with the current
    % PREFS, as the discretization might have changed in the periodic case:
    varargin{nargin-1} = bcType;
    varargin{nargin} = prefs;
else
    % Otherwise, add the PREFS to VARARGIN, so that it can be passed to the call
    % to LINOP/SVDS below.
    varargin{nargin} = bcType;
    varargin{nargin+1} = prefs;
end

% construct adjoint
[varargout{1:nargout}] = svds(L,varargin{:});

end