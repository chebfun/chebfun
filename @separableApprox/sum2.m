function I = sum2( f, varargin )
%SUM2   Double integral of a SEPARABLEAPPROX over its domain. 
%   I = SUM2(F) returns the double integral of a SEPARABLEAPPROX. 
% 
% See also INTEGRAL2, INTEGRAL, QUAD2D. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

I = integral2( f, varargin{:} ); 

end
