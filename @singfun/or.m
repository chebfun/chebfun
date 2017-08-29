function h = or(f, g)
%|   SINGFUN logical OR.
%   F | G performs a logical OR of the SINGFUN objects F and G and returns a
%   SMOOTHFUN containing elements set to either logical 1 (TRUE) or logical 0
%   (FALSE).  The output SMOOTHFUN is set to 1 if either input SINGFUN
%   contains a non-zero element at that point, otherwise it is set to 0.
%   F and G must either be identically zero or have roots in their domains.
%   If this is not the case, garbage is returned with no warning.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'singfun') )
    f = f.smoothPart;
end

if ( isa(g, 'singfun') )
    g = g.smoothPart;
end
    
h = f | g;

end
