function I = integral( f, varargin )
%INTEGRAL Complete definite integral of chebfun2. 
%
% I = INTEGRAL(F), returns the definite integral of a chebfun2. Integrated
% over its domain of definition.
% 
% I = INTEGRAL(F,g), returns the integral of a chebfun2 along the curve
% defined by the complex-valued chebfun g. 
% 
% See also INTEGRAL2, SUM2, QUAD2D.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin == 1 )                        % another way to do sum2(f) 
   I = integral2( f ); 
else
    if ( ~isa( varargin{1}, 'chebfun' ) )  % Line integral over a chebfun
        I = integral2( f, varargin{ : } );
    else                                  
        % Get curve: 
        c = varargin{1}; 
        
        % make complex: 
        c = c + realmin*1i;   
        
        % line integral: 
        I = sum( feval(f, c ) .* abs( diff( c ) ) );
    end
end

end