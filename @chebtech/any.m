function a = any(f, dim)
%ANY   True if any element of a CHEBTECH is a nonzero number. ANY ignores
%      entries that are NaN (Not a Number).
%   ANY(F, DIM), where F is an array-valued CHEBTECH, works down the dimension
%   DIM.  If DIM is 1, then ANY returns a logical row vector in which the Jth
%   element is TRUE if any element of the Jth column is nonzero.  If DIM is 2,
%   ANY returns a CHEBTECH which takes the value 1 wherever any of the columns
%   (or rows) of F are nonzero, and zero everywhere else.  In this case, F must
%   either be identically zero or have no roots in its domain.  Otherwise,
%   garbage is returned without warning.
%
%   ANY(F) is shorthand for ANY(F, 1).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information

% Parse inputs:
if ( nargin < 2 )
    dim = 1;
end

if ( dim  == 1 )        % ANY down the columns.
    a = any(f.coeffs);
elseif ( dim == 2 )     % ANY down the rows.
    a = f;
    arbitraryPoint = 0.1273881594;
    a.coeffs = any(feval(a, arbitraryPoint));
else
    error('CHEBFUN:CHEBTECH:any:dim', 'DIM input must be 1 or 2.');
end

end
