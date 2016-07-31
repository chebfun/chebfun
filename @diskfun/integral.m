function I = integral( f, varargin )
%INTEGRAL   Complete definite integral of DISKFUN. 
%
%   I = INTEGRAL(F), returns the definite integral of a DISKFUN integrated
%   over its domain of definition.
% 
%   I = INTEGRAL(F, g), returns the integral of F along the 
%   curve defined by the complex-valued CHEBFUN g.
% 
%   I = INTEGRAL(F, 'unitcircle') returns the integral of F along the 
%   unit circle. 
% See also SUM2. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )                         % Another way to do sum2(f) 
    
   I = integral2( f ); 
   
else  
    if strcmpi(varargin{1},'unitcircle')
            c = chebfun(@(t) exp(1i*t), [-pi, pi]); 
    else
             c = varargin{1}; 
    end
    if ( ~isa( c, 'chebfun' ) )  % Line integral over a CHEBFUN
        I = integral2( f, varargin{ : } );
    else                                  
        % Get curve: 
        
        %c = varargin{1}; 
        
        % Make complex: 
        c = c + realmin*1i;   
        % Line integral: 
        I = sum( feval(f, c ) .* abs( diff( c ) ) );
    end
    
    
end

end