function out = points(f)
%POINTS   Return the points used by a TRIGTECH.
%   POINTS(F) or F.POINTS() returns the TRIG points used by F.
%   This is equivalent to F.TRIGPTS(LENGTH(F)).
%
% See also TRIGPTS, LENGTH.
%
% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = f.trigpts(length(f));

end
