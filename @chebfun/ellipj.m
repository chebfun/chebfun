function [sn, cn, dn] = ellipj(u, m, pref)
%ELLIPJ Jacobi elliptic functions.
%   [SN, CN, DN] = ELLIPJ(U, M) returns the chebfuns of the Jacobi elliptic
%   functions Sn, Cn, and Dn with parameter M composeosed with the chebfun U.M must
%   be a scalar in the range 0 <= M <= 1.
%
%   [SN, CN, DN] = ELLIPJ(U, M, TOL) composeutes the elliptic functions to the
%   accuracy TOL instead of the default TOL = CHEBFUN.PREF('EPS').
%
%   Complex values of U are accepted, but the resulting composeutation may be
%   inaccurate. Use ELLIPJC from Driscoll's SC toolbox instead.
%
%   Note:Some definitions of the Jacobi elliptic functions use the modulus k
%   instead of the parameter M. They are related by M = k^2.
%
% See also ELLIPKE

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin < 3 )
    pref = chebfun.pref();
    tol = pref.chebfun.eps;
elseif ( ~isstruct(pref) )
    tol = pref;
    pref = chebfun.pref();
else
    tol = pref.chebfun.eps;
end

if ( isreal(u) )

    % SN
    sn = compose(u, @(x) ellipj(x, m, tol), pref);

    % CN
    if ( nargout >= 2 )
        cn = compose(u, @(x) cnfun(x, m, tol), pref);
    end

    % DN
    if ( nargout == 3 )
        dn = compose(u, @(x) dnfun(x, m, tol), pref);
    end

else
    % Use imaginary transformations:

    [s, c, d] = ellipj(real(u), m, tol);      % real values
    [s1, c1, d1] = ellipj(imag(u), 1-m, tol); % imaginary values
    denom = c1.^2 + m*(s.*s1).^2;

    % SN
    sn = (s.*d1 + 1i*c.*d.*s1.*c1)./denom;

    % CN
    if ( nargout >= 2 )
        cn = (c.*c1-1i*s.*d.*s1.*d1)./denom;
    end

    % DN
    if ( nargout == 3 )
        dn = (d.*c1.*d1-1i*m*s.*c.*s1)./denom;
    end

end

end

function cnout = cnfun(u, varargin)
    [ignored, cnout, ~] = ellipj(u,varargin{:});
end

function dnout = dnfun(u, varargin)
    [ignored, ignored, dnout] = ellipj(u, varargin{:});
end
