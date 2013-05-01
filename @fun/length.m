function len = length(f)
%LENGTH	  Length of a BNDFUN.
%   LENGTH(F) is the number of values at Chebyshev points used to represent F.
%
%   See also SIZE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The length of a BNDFUN is the length of its onefun.
len = length(f.onefun);

end
