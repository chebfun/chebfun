function g = cos( f )
%COS   Cosine of a BALLFUN.
%   COS(F) computes the cosine of F.
%
% See also COSH, COS.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = compose( f, @cos ); 

end
