function N = normal( F, varargin )
%NORMAL normal vector to a surface represented by a chebfun2v.
%   N = NORMAL(F) returns a chebfun2v representing the normal vector to the
%   surface F. The vector has the same magntiude as the surface's tangent vector
%
%   N = NORMAL(F,'unit') returns the unit normal vector, represented as a
%   chebfun2v, to the surface F.
%
% See also CHEBFUN/NORMAL. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

N = cross( diff( F, 1, 2 ), diff( F, 1, 1 ) );

if ( nargin > 1 )
    if ( strcmpi(varargin{1}, 'unit') )
        if ( N.nComponents == 2 )
            r = roots( N );
            if ( ( ~isempty(r) ) || ( norm(N) == 0 ) )
                error('CHEBFUN:NORMAL:ZERO', 'Normal vector is zero.');
            end
            N = N ./ abs( N );
        elseif ( norm(N) == 0 )
            error('CHEBFUN:NORMAL:ZERO', 'Normal vector is zero');
        end          
        
    else
        error('CHEBFUN:NORMAL', 'Second argument is not recognised.');
    end
end

end