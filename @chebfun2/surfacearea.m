function I = surfacearea(f , varargin )
%SURFACEAREA    Surface area of a CHEBFUN2.
%   SURFACEAREA(F) computes the surface area of the CHEBFUN2 in the domain of F.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin > 1 )
    if ( isa(varargin{1}, 'double') )
        if ( length(varargin{1}) == 4 )
            % Restricted surface area over rectangular region. 
            dom = varargin{1}; 
            f = restrict( f, dom );
        else
            error('CHEBFUN:CHEBFUN2:surfacearea:domain', 'Bad domain.');
        end
    elseif ( isa(varargin{1}, 'chebfun') )
        f = restrict( f, varargin{1} );
        % Surface area is now just the arc length. 
        I = sum( sqrt( 1 + diff( f ).^2 ) ); 
        return
    else
        error('CHEBFUN:CHEBFUN2:surfacearea:domain','Bad restricting domain.');
    end
end

% First order derivatives:
fx = diff(f, 1, 2); 
fy = diff(f, 1, 1);  

% Integrand:
G = 1 + fx.^2 + fy.^2;
S = sqrt( G );

% Surface area:
I = integral2( S );

end
