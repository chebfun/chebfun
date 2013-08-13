function out = length(f)
%LENGTH   Length of a Chebfun.
%   LENGTH(F) returns the length of a CHEBFUN object F, which is defined as the
%   sum of the length of f.funs.
%
% See also CHEBFUN/SIZE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = sum(cellfun(@length, f.funs));

end