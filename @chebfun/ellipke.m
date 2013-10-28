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
%   See also ELLIPJ.

% Choose a tolerance:
tol = get(m, 'epslevel');
if ( nargin == 1 )
    pref = chebpref();
    tol = max(pref.techPrefs.eps, tol);
elseif ( isnumeric(pref) )
    tol = max(pref, tol);
    pref = chebpref();
    pref.techPrefs.eps = tol;
else
    tol = max(pref.techPrefs.eps, tol);
end

% Call COMPOSE():
try 
    k = compose(m, @(m) ellipke(m, .1*tol), pref);
catch ME
   if ( strcmp(ME.identifier, 'MATLAB:ellipke:MOutOfRange') )
       error('CHEBFUN:ellipke:MOutOfRange', ...
           'M must be in the range 0 <= M <= 1.');
   else
       rethrow(ME)
   end
end
    
% Compute the second complete elliptic integral if required:
if ( nargout == 2 )
    % Call COMPOSE():
    e = compose(m, @(m) eFun(m, tol), pref);
end

end

function e = eFun(m, tol)
    [ignored, e] = ellipke(m, tol);
end
