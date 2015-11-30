function [sn, cn, dn] = ellipj(u, m, pref)
%ELLIPJ   Jacobi elliptic functions.
%   [SN, CN, DN] = ELLIPJ(U, M) returns CHEBFUNS for the compositions Sn(U)
%   Cn(U), and Dn(U), where Sn, Cn, and Dn are the Jacobi elliptic functions
%   with parameter M. U may be a scalar or a CHEBFUN, and M must be a CHEBFUN
%   or scalar in the range 0 <= M <= 1.
%
%   [SN, CN, DN] = ELLIPJ(U, M, TOL) composes the elliptic functions to the
%   accuracy TOL instead of the default TOL = EPS.
%
%   Complex values of U are accepted, but the resulting computation may be
%   inaccurate. Use ELLIPJC from Driscoll's SC toolbox instead.
%
%   Note: Some definitions of the Jacobi elliptic functions use the modulus k
%   instead of the parameter M. They are related by M = k^2.
%
% See also ELLIPKE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note. It isn't possible to call MATLAB's built in cn and/or dn in isolation.
% To get around these we use the subfunctions cnFun and dnFun.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( nargin < 3 )
    pref = chebfunpref();
    tol = pref.techPrefs.eps;
elseif ( ~isstruct(pref) )
    tol = pref;
    pref = chebfunpref();
else
    tol = pref.techPrefs.eps;
end

    function x = fudge(x, tol)
        % M must lie in [0, 1]. Fudge values outside this range during compose.
        x(x < 0 & x > -tol,:) = 0;
        x(x > 1 & x < 1 + tol,:) = 0;
    end

if ( isreal(u) )
    % Real case: Call ELLIPJ().

    if ( isnumeric(m) )     % U = CHEBFUN, M = double
        % SN
        sn = compose(u, @(u) ellipj(u, m, tol), pref);
        % CN
        if ( nargout >= 2 )
            cn = compose(u, @(u) cnFun(u, m, tol), pref);
        end
        % DN
        if ( nargout == 3 )
            dn = compose(u, @(u) dnFun(u, m, tol), pref);
        end
    elseif ( isnumeric(u) )  % U = double, M = CHEBFUN
        mTol = max(eps*vscale(m), tol);
        % SN
        sn = compose(m, @(m) ellipj(u, fudge(m, mTol), tol), pref);
        % CN
        if ( nargout >= 2 )
            cn = compose(m, @(m) cnFun(u, fudge(m, mTol), tol), pref);
        end
        % DN
        if ( nargout == 3 )
            dn = compose(m, @(m) dnFun(u, fudge(m, mTol), tol), pref);
        end
    else                     % U = CHEBFUN, M = CHEBFUN
        mTol = max(eps*vscale(m), tol);
        % SN
        sn = compose(u, @(u, m) ellipj(u, fudge(m, mTol), tol), m, pref);
        % CN
        if ( nargout >= 2 )
            cn = compose(u, @(u, m) cnFun(u, fudge(m, mTol), tol), m, pref);
        end
        % DN
        if ( nargout == 3 )
            dn = compose(u, @(u, m) dnFun(u, fudge(m, mTol), tol), m, pref);
        end
    end
    
else
    % Use imaginary transformations (see http://dlmf.nist.gov/22.8):

    [s, c, d] = ellipj(real(u), m, tol);      % real values
    [s1, c1, d1] = ellipj(imag(u), 1-m, tol); % imaginary values
    denom = 1./(c1.^2 + m.*(s.*s1).^2);

    % SN ( see NIST 22.8.1)
    sn = (s.*d1 + 1i*c.*d.*s1.*c1).*denom;

    % CN ( see NIST 22.8.2)
    if ( nargout >= 2 )
        cn = (c.*c1 - 1i*s.*d.*s1.*d1).*denom;
    end

    % DN ( see NIST 22.8.3)
    if ( nargout == 3 )
        dn = (d.*c1.*d1 - 1i*m.*s.*c.*s1).*denom;
    end

end

end

function cnOut = cnFun(u, varargin)
%CNFUN  Jacobi elliptic function cn.
    [ignored, cnOut, ignored] = ellipj(u, varargin{:});
end

function dnOut = dnFun(u, varargin)
%DNFUN  Jacobi elliptic function dn.
    [ignored, ignored, dnOut] = ellipj(u, varargin{:});
end


