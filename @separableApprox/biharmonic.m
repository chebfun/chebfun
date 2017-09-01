function B = biharmonic(f)
%BIHARMONIC   Biharmonic operator of a SEPARABLEAPPROX.
%   B = BIHARMONIC(F) returns a SEPARABLEAPPROX representing the biharmonic 
%   operator applied to F.
%
% See also SEPARABLEAPPROX/BIHARM.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% biharmonic(f) = f_xxxx + f_yyyy + 2*f_xxyy:
B = diff(f, 4, 2) + diff(f, 4, 1) + 2*diff(diff(f, 2, 1), 2, 2);

end
