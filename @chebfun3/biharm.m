function B = biharm(f)
%BIHARM   Biharmonic operator applied to a CHEBFUN3.
%   B = BIHARM(F) returns a CHEBFUN3 object B representing the biharmonic 
%   operator applied to a CHEBFUN3 object F.
%
%   This is shorthand for BIHARMONIC(F).
%
% See also CHEBFUN3/BIHARMONIC.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call Biharmonic: 
B = biharmonic(f);

end