function out = points(f)
%POINTS   Return the points used by a CHEBTECH.
%   POINTS(F) or F.POINTS() returns the Chebyshev points used by F.
%   This is equivalent to F.CHEBPTS(LENGTH(F)).
%
% See also CHEBPTS, LENGTH.
%
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = f.chebpts(length(f));

end
