function len = length(f)
%LENGTH	  Length of a CLASSICFUN.
%   LENGTH(F), where F is a CLASSICFUN, is the length of the ONEFUN of F.
%
% See also SIZE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The length of a CLASSICFUN is the length of its ONEFUN.
len = length(f.onefun);

end
