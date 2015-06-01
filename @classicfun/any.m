function a = any(f, dim)
%ANY   True if any element of a CLASSICFUN is a nonzero number. ANY ignores
%      entries that are NaN (Not a Number).
%   ANY(X, DIM), where X is an array-valued CLASSICFUN, works down the dimension DIM.
%   If DIM is 1, then ANY returns a logical row vector in which the Jth element
%   is TRUE if any element of the Jth column is nonzero.  If DIM is 2, ANY
%   returns a CLASSICFUN which takes the value 1 wherever any of the columns (or rows)
%   of X are nonzero, and zero everywhere else.  In this case, X must either be
%   identically zero or have no roots in its domain.  Otherwise, garbage is
%   returned without warning.
%
%   ANY(X) is shorthand for ANY(X, 1).

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information

if ( nargin < 2 )
    dim = 1;
end

if ( dim == 1 )
    a = any(f.onefun, 1);
elseif ( dim == 2 )
    a = f;
    a.onefun = any(a.onefun, 2);
else
    error('CHEBFUN:CLASSICFUN:any:dim', 'DIM input must be 1 or 2.');
end

end
