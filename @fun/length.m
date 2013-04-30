function len = length(f)
%LENGTH	  Length of a FUN.
%   LENGTH(F), where F is a FUN, is the length of the ONEFUN of F.
%
%   See also SIZE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The length of a FUN is the length of its onefun.
len = length(f.onefun);

end
