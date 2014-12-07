function [f, g] = transferDeltas(f, g)
%TRANSFERDELTAS   Transfers delta functions in F that are at the right end
%   F to G.
%
% See also DELTAFUN/MERGEDELTAS

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If F doesn't have delta functions, return:
if ( isempty(f.deltaLoc) )
    return;
end

% Get the tolerance:
pref = chebfunpref();
tol = pref.deltaPrefs.proximityTol;

% Get the domains:
dom = f.domain;
b = dom(end);
dom = g.domain;
c = dom(1);

% If the domains are not contiguous, return:
if ( abs(c - b) > tol )
    return
end

% If F doesn't have a delta function at the right end point, return:
if ( abs(f.deltaLoc(end) - c) > tol )
    return
end

% F does have a delta function at the right end point, transfer it to G:

% Copy the delta function and remove it from F:
deltaRight = f.deltaMag(:, end);
f.deltaMag(:, end) = [];
f.deltaLoc(end) = [];

% Transfer it to G:
if ( isa(g, 'deltafun') )
    [g.deltaMag, g.deltaLoc] = deltafun.mergeDeltas(deltaRight, c, g.deltaMag, g.deltaLoc);
else
    g = deltafun(g, struct('deltaMag', deltaRight, 'deltaLoc', c));
end
f = simplifyDeltas(f);
    
