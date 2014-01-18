function I = sum2( f, varargin )
%SUM2 Double integral of a chebfun2 over its domain. 
%
% I = SUM2(F) returns the double integral of a chebfun2. 
% 
% See also INTEGRAL2, INTEGRAL, QUAD2D. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

I = integral2( f, varargin{:} ); 

end
