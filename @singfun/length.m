function len = length(f)
%LENGTH   Length of a SINGFUN.
%   LENGTH(F) is the length of the smooth part contained in F.
%
% See also SIZE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The length of a SINGFUN object is the length of its smooth part:
len = length(f.smoothPart);

end
