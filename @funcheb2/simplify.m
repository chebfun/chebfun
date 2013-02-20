function f = simplify(f, pref)
%SIMPLIFY  Trim trailing Chebyshev coefficients of a funcheb2. 
%   G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the
%   funcheb2 F such that length(G) <= length(F) but ||G - F|| < F.epslevel.
%
%   G = SIMPLIFY(F, PREF) does the same, but using some of the preferences in
%   the preference structure PREF. In particular, by default SIMPLIFY will use
%   funcheb2.classicCheck but an alternative can be passed in the
%   PREF.funcheb2.happinessCheck field. Additionally, G will be computed to
%   satisfy ||G - F|| < max(F.epslevel, PREF.funcheb2.eps).
%
% See also funcheb2.happinessCheck.m, funcheb2.pref.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Grab some preferences:
if ( nargin < 2 )
    pref = funcheb2.pref;
end

% Take max of pref.eps and epslevel:
pref.funcheb2.eps = max(pref.funcheb2.eps, f.epslevel);

% Check to see if we can trim the tail:
[ishappy, epslevel, cutoff] = ...
    funcheb2.happinessCheck([], f.values, f.coeffs, f.vscale, pref);

% Trim/alias it with prolong:
f.ishappy = ishappy | f.ishappy;
f.epslevel = epslevel;
f = prolong(f, cutoff); 
    
end
