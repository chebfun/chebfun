function g = sign(f)
%SIGN    Sign function of a CHEBFUN.
%   G = SIGN(F) returns a piecewise constant CHEBFUN G such that G(x) = 1 in the
%   interval where F(x) > 0, G(x) = -1 in the interval where F(x) < 0 and G(x) =
%   0 in the interval where F(x) = 0. Breakpoints in G are introduced at zeros
%   of F.
%
%   For the nonzero elements of complex F, SIGN(F) = F./ABS(F)
%
% See also ABS, HEAVISIDE, ROOTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Deal with the empty case:
if ( isempty(f) )
    g = f;
    return
end

% Add breaks at the appropriate roots of f:
g = addBreaksAtRoots(f);

% Call SIGN on each of the FUNs: (result will be smooth)
for k = 1:numel(g.funs)
    g.funs{k} = sign(g.funs{k});
end

% Take the absolute value of the impulses in the first row:
g.impulses = sign(g.impulses(:,:,1));

% Remove unnecessary break points:
g = merge(g); 

end
