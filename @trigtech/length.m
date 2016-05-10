function len = length(f)
%LENGTH   Length of a TRIGTECH.
%   LENGTH(F) is the number of values at equally spaced used to represent F.
%
% See also SIZE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The length of a TRIGTECH is the length of its vector of values.
len = size(f.values, 1);

end
