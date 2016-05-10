function f = fliplr(f)
%FLIPLR    Flip columns of an array-valued CLASSICFUN object.
%   FLIPLR(F) flips the columns of an array-valued CLASSICFUN F in the left/right
%   direction. If F has only one column, then this function has no effect.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Flip the ONEFUN:
f.onefun = fliplr(f.onefun);

end
