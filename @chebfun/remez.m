function varargout = remez(f, varargin)
%REMEZ   Best polynomial or rational approximation for real valued chebfuns.
%   P = REMEZ(F, M) computes the minimax polynomial approximation of degree M
%   to the real CHEBFUN F using the Remez algorithm.
%
%   [P, Q] = REMEZ(F, M, N) computes the minimax rational approximation P/Q
%   of type (M, N).
%
%   [P, Q, R_HANDLE] = REMEZ(F, M, N) additionally returns a function handle
%   for evaluating P/Q.
%
%   [...] = REMEZ(..., 'tol', TOL) uses the value TOL as the termination
%   tolerance on the relative equioscillation error.
%
%   [...] = REMEZ(..., 'display', 'iter') displays output at each iteration.
%
%   [...] = REMEZ(..., 'maxiter', MAXITER) sets the maximum number of allowable
%   iterations to MAXITER.
%
%   [...] = REMEZ(..., 'init', XK) allows the user to specify the vector XK as
%   the starting reference.
%
%   [...] = REMEZ(..., 'plotfcns', 'error') plots the error after each iteration
%   while the algorithm executes.
%
%   [P, ERR] = REMEZ(...) and [P, Q, R_HANDLE, ERR] = REMEZ(...) returns the
%   maximum error ERR.
%
%   [P, ERR, STATUS] = REMEZ(...) and [P, Q, R_HANDLE, ERR, STATUS] = REMEZ(...)
%   return a structure array STATUS with the following fields:
%      STATUS.DELTA  - Obtained tolerance.
%      STATUS.ITER   - Number of iterations performed.
%      STATUS.DIFFX  - Maximum correction in last trial reference.
%      STATUS.XK     - Last trial reference on which the error equioscillates.
%
%   This code is quite reliable for polynomial approximations but may sometimes
%   have difficulties in the rational case.
%
% References:
%
%   [1] B. Beckermann, S. Filip and Y. Nakatsukasa, manuscript in preparation.
%
%   [2] R. Pachon and L. N. Trefethen, "Barycentric-Remez algorithms for best
%   polynomial approximation in the chebfun system", BIT Numerical Mathematics,
%   49:721-742, 2009.
%
%   [3] R. Pachon, "Algorithms for Polynomial and Rational Approximation".
%   D. Phil. Thesis, University of Oxford, 2010 (Chapter 6).
%
% See also CF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:CHEBFUN:remez', ...
        [' This command is deprecated.', ...
         ' The user should use minimax instead.']);
polyOutput = detectType(f,varargin{:});

if ( polyOutput )
    [p,err,status] = minimax(f,varargin{:});
    varargout = {p, err, status};
else
    [p,q,r,err,status] = minimax(f,varargin{:});
    varargout = {p, q, r, err, status};
end

end
    
function polyOutput = detectType(~, varargin)

isSilent = 0;
for k = 1:length(varargin)
    if ( ischar(varargin{k}) && strcmpi('silent', varargin{k}) )
        isSilent = 1;
    end
end

% Detect polynomial / rational approximation type.
polyOutput = true;
if ( mod(nargin - isSilent, 2) ) % Odd number of inputs --> rational case.                             
    polyOutput = false;
end

end