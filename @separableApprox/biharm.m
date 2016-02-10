function B = biharm(f)
%BIHARM   Biharmonic operator of a SEPARABLEAPPROX.
%   B = BIHARM(F) returns a SEPARABLEAPPROX representing the biharmonic operator 
%   applied to F.
%
%   This is shorthand for BIHARMONIC(F).
%
% See also SEPARABLEAPPROX/BIHARMONIC.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call SEPARABLEAPPROX/BIHARMOMIC:
B = biharmonic(f);

end
