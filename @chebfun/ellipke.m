function [k, e] = ellipke(m, pref)
%ELLIPKE   Complete elliptic integral of a CHEBFUN.
%   [K, E] = ELLIPKE(M) returns the value of the complete elliptic integrals of
%   the first and second kinds, composed with the CHEBFUN M.  As currently
%   implemented, M is limited to 0 <= M <= 1.
%
%   [K, E] = ELLIPKE(M, TOL) computes the complete elliptic integrals to the
%   accuracy TOL instead of the default TOL = EPS(CLASS(M)).
%
%   Some definitions of the complete elliptic integrals use the modulus k
%   instead of the parameter M.  They are related by M = k^2.
%
% See also ELLIPJ.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Choose a tolerance:
tol = eps;
if ( nargin == 1 )
    pref = chebfunpref();
    tol = max(pref.techPrefs.eps, tol);
elseif ( isnumeric(pref) )
    tol = max(pref, tol);
    pref = chebfunpref();
    pref.techPrefs.eps = tol;
else
    tol = max(pref.techPrefs.eps, tol);
end

    function x = fudge(x, tol)
        % M must lie in [0, 1]. Fudge values outside this range during compose.
        x(x < 0 & x > -tol,:) = 0;
        x(x > 1 & x < 1 + tol,:) = 0;
    end

% Loop over the columns:
for j = numel(m):-1:1
    mTol = max(eps*vscale(m(j)), tol);
    try
        % Call COMPOSE():
        k(j) = compose(m(j), @(m) ellipke(fudge(m, mTol), .1*tol), pref);
    catch ME
        if ( strcmp(ME.identifier, 'MATLAB:ellipke:MOutOfRange') )
            error('CHEBFUN:CHEBFUN:ellipke:MOutOfRange', ...
                'M must be in the range 0 <= M <= 1.');
        else
            rethrow(ME)
        end
    end
end

% Compute the second complete elliptic integral if required:
if ( nargout == 2 )
    e(numel(m)) = chebfun();
    for j = numel(m):-1:1
        mTol = max(eps*vscale(m(j)), tol);
        % Call COMPOSE():
        e(j) = compose(m(j), @(m) eFun(fudge(m, mTol), tol), pref);
    end
end

end

function e = eFun(m, tol)
    [ignored, e] = ellipke(m, tol);
end
