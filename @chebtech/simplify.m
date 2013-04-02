function f = simplify(f, pref, type)
%SIMPLIFY   Trim trailing Chebyshev coefficients of a CHEBTECH object. 
%   G = SIMPLIFY(F) attempts to compute a 'simplified' version G of the CHEBTECH
%   object F such that LENGTH(G) <= LENGTH(F) but ||G - F|| < F.EPSLEVEL.
%
%   G = SIMPLIFY(F, PREF) does the same, but using some of the preferences in
%   the preference structure PREF. In particular, SIMPLIFY will use CLASSICCHECK
%   to test for accuracy by default, but an alternative can be passed in the
%   PREF.CHEBTECH.HAPPINESSCHECK field. Additionally, G will be computed to
%   satisfy ||G - F|| < MAX(F.EPSLEVEL, PREF.CHEBTECH.EPS).
%
% See also HAPPINESSCHECK, CLASSICCHECK, PREF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Deal with empty case:
if ( isempty(f) )
    return
end

% Grab some preferences:
if ( nargin < 2 )
    pref = chebtech.pref();
elseif ( isnumeric(pref) )
    pref = chebtech.pref('eps', pref);
elseif ( ~isstruct(pref) )
    pref = chebtech.pref('happinessCheck', pref);
end
if ( nargin == 3 )
    pref.chebtech.happinessCheck = type;
end

% Take max of PREF.EPS and EPSLEVEL:
pref.chebtech.eps = max(pref.chebtech.eps, f.epslevel);

% Check to see if we can trim the tail:
[ishappy, epslevel, cutoff] = happinessCheck(f, pref);
% Trim/alias it with prolong:
f.ishappy = ishappy || f.ishappy;
f.epslevel = epslevel;
f = prolong(f, cutoff);
 
end
