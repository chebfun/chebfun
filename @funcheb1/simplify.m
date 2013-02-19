function f = simplify(f, pref)
%SIMPLIFY  Trim trailing Chebyshev coefficients of a funcheb1. 
%   G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the
%   funcheb1 F such that length(G) <= length(F) but ||G - F|| < F.epslevel.
%
%   G = SIMPLIFY(F, PREF) does the same, but using some of the preferences in
%   the preference structure PREF. In particular, by default SIMPLIFY will use
%   funcheb1.checkCoeffs but an alternative can be passed in the
%   PREF.funcheb1.happinessCheck field. Additionally, G will be computed to
%   satisfy ||G - F|| < max(F.epslevel, PREF.funcheb1.eps).
%
% See also funcheb1.isHappy.m, funcheb1.pref.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Grab some preferences:
if ( nargin < 2 )
    pref = funcheb1.pref;
end

% Take max of pref.eps and epslevel:
pref.funcheb1.eps = max(pref.funcheb1.eps, f.epslevel);

% Check to see if we can trim the tail:
if ( strcmpi(pref.funcheb1.happinessCheck, 'default') )
    % Call default happiness check (but without resampling):
    [cutoff, epslevel] = funcheb1.checkCoeffs(f.values, f.coeffs, f.vscale, pref);
else
    [cutoff, epslevel] = pref.funcheb1.happinessCheck(@(x) feval(f, x), ...
        f.values, f.coeffs, f.vscale, pref);
end

% Trim/alias it with prolong:
f.epslevel = epslevel;
f = prolong(f, cutoff); 
    
end
