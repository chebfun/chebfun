function out = points(f)
%POINTS   Return the points used by a FOURTECH.
%   POINTS(F) or F.POINTS() returns the Fourier points used by F.
%   This is equivalent to F.FOURPTS(LENGTH(F)).
%
% See also FOURPTS, LENGTH.
%
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = f.fourpts(length(f));

end