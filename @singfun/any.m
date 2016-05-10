function a = any(f, dim)
%ANY   True if any element of a SINGFUN is a nonzero number. ANY ignores
%      entries that are NaN (Not a Number).
%   ANY(F, DIM), where DIM is the dimension along which ANY is taken. Note
%   that array-valued SINGFUN is not supported. If DIM is 1, then ANY 
%   returns a logical value which is TRUE if any element of F is nonzero.  
%   If DIM is 2, ANY returns a SMOOTHFUN which takes the value 1 if F is
%   nonzero.  In this case, F must either be identically zero or have no 
%   roots in its domain.  Otherwise, garbage is returned without warning.
%
%   ANY(F) is shorthand for ANY(F, 1).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information

% Parse inputs:
if ( nargin < 2 )
    dim = 1;
end

a = any(f.smoothPart, dim);

end
