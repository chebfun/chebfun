function h = and(f, g)
%&   CHEBFUN Logical AND
%   F & G performs a logical AND of two CHEBFUN objects F and G and returns a
%   CHEBFUN with values set to either logical 1 (TRUE) or logical 0 (FALSE).  A
%   value of the output CHEBFUN is set to 0 if both input CHEBFUN objects have
%   a non-zero value at that point, otherwise it is set to 0.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with the empty case:
if ( isempty(f) || isempty(g) )
    h = chebfun();
    return
end

% Make sure dimensions and domains are the same:
if ( ~isequal(size(f), size(g)) )
    error('CHEBFUN:CHEBFUN:and:dims', 'Dimensions of f and g must agree.');
elseif ( ~domainCheck(f, g) )
    error('CHEBFUN:CHEBFUN:and:doms', 'Domains of f and g must agree.');
end

% Split up f and g at their roots, and ensure breakpoints are the same:
f = addBreaksAtRoots(f);
g = addBreaksAtRoots(g);
[f, g] = overlap(f, g);

% Call AND() on the FUNs:
h = f;
for k = 1:1:numel(h.funs)
    h.funs{k} = f.funs{k} & g.funs{k};
end

% Deal with point values:
h.pointValues = f.pointValues & g.pointValues;

% Get rid of unnecessary breakpoints:
h = merge(h);

end
