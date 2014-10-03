function len = length(f)
%LENGTH   Length of a FOURTECH.
%   LENGTH(F) is the number of values at Fourier points (equally spaced)
%   used to represent F.
%
% See also SIZE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The length of a FOURTECH is the length of its vector of values.
len = size(f.values, 1);

end