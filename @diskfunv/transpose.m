function F = transpose( F )
% .' transpose of a DISKFUN
%    F = TRANSPOSE(F) is the same as F.'. The orientation of the DISKFUNV 
%    is transposed. 
%
%    See also DISKFUNV/CTRANSPOSE 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

F.isTransposed = ~F.isTransposed;

end