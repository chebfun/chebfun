function X = mldivide(A, B)
%\   Left matrix divide for SINGFUN objects.
%
%   Since SINGFUN does not support array-valued objects, left matrix divide is 
%   exactly same as left divide.
%
% See also LDIVIDE, RDIVIDE, MRDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

X = ldivide(A, B);

end
