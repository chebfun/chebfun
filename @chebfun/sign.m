function F = sign(F, pref)
%SIGN   Sign function of a CHEBFUN.
%   G = SIGN(F) returns a piecewise constant CHEBFUN G such that G(x) = 1 in the
%   interval where F(x) > 0, G(x) = -1 in the interval where F(x) < 0 and G(x) =
%   0 in the interval where F(x) = 0. Breakpoints in G are introduced at zeros
%   of F.
%
%   For the nonzero values of complex F, SIGN(F) = F./ABS(F)
%
% See also ABS, HEAVISIDE, ROOTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with the empty case:
if ( isempty(F) )
    return
end

if ( nargin < 2 ) 
    pref = chebfunpref();
end

for k = 1:numel(F)
    F(k) = signColumn(F(k), pref);
end

end

function g = signColumn(f, pref)

% Add breaks at the appropriate roots of f:
g = addBreaksAtRoots(f, pref);

% Call SIGN on each of the FUNs: (result will be smooth)
for k = 1:numel(g.funs)
    g.funs{k} = sign(g.funs{k}, pref);
end

% Take the sign of the pointValues in the first row:
g.pointValues = sign(g.pointValues);

% Remove unnecessary breakpoints:
g = merge(g, pref); 

end
