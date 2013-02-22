function f = simplify(f, pref)
%SIMPLIFY  Trim trailing Chebyshev coefficients of a FUNCHEB object. 
%   G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the FUNCHEB
%   object F such that length(G) <= length(F) but ||G - F|| < F.epslevel.
%
%   G = SIMPLIFY(F, PREF) does the same, but using some of the preferences in
%   the preference structure PREF. In particular, SIMPLIFY will use classicCheck
%   by default, but an alternative can be passed in the
%   PREF.(class(f)).happinessCheck field. Additionally, G will be computed to
%   satisfy ||G - F|| < max(F.epslevel, PREF.(class(f)).eps).
%
% See also funcheb.happinessCheck, classicCehck, funcheb.pref.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Grab some preferences:
if ( nargin < 2 )
    pref = f.pref();
end

% Take max of pref.eps and epslevel:
pref.(class(f)).eps = max(pref.(class(f)).eps, f.epslevel);

% Check to see if we can trim the tail:
[ishappy, epslevel, cutoff] = happinessCheck(f, pref);

% Trim/alias it with prolong:
f.ishappy = ishappy | f.ishappy;
f.epslevel = epslevel;
f = prolong(f, cutoff); 
    
end
