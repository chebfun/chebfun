function biharmF = biharmonic(f)
%BIHARMONIC   Biharmonic operator of a DISKFUN.
%   B = BIHARMONIC(F) returns a DISKFUN representing the biharmonic 
%   operator applied to F.
%
% See also DISKFUN/BIHARM.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Use the implementation in biharm.
biharmF = biharm(f);

end
