function biharmF = biharm(f)
%BIHARM   Biharmonic operator of a SPHEREFUN.
%   B = BIHARM(F) returns a SPHEREFUN representing the biharmonic operator 
%   applied to F.
%
%   This is shorthand for BIHARMONIC(F).
%
% See also SPHEREFUN/BIHARMONIC.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Implement as the laplacian of the laplacian
biharmF = laplacian(laplacian(f));

end
