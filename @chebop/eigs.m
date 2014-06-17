function varargout = eigs(N, varargin)
%EIGS   Find selected eigenvalues and eigenfunctions of a linear CHEBOP.
%   D = EIGS(A) returns a vector of 6 eigenvalues of the linear CHEBOP A. EIGS
%   will attempt to return the eigenvalues corresponding to the least
%   oscillatory eigenfunctions. (This is unlike the built-in EIGS, which returns
%   the largest eigenvalues by default.). If A is not linear, an error is
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

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Linearize and check whether the chebop is linear:
[L, ignored, fail] = linop(N); %#ok<ASGLU>

% Support for generalised problems:
if ( ~fail && nargin > 1 && isa(varargin{1}, 'chebop') )
    % Linearise the second CHEBOP:
    [varargin{1}, ignored, fail] = linop(varargin{1}); %#ok<ASGLU>
end

if ( fail )
    error('CHEBFUN:CHEBOP:eigs:nonlinear', ...
        ['The operator appears to be nonlinear.\n', ...
         'EIGS() supports only linear CHEBOP instances.']);
end

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
