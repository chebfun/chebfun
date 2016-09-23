function len = length(f)
%LENGTH   Length of a DELTAFUN.
%   LENGTH(F) is the length of the smooth part contained in F.
%
% See also SIZE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The length of a DELTAFUN object is the length of its funPart.
len = length(f.funPart);

end
