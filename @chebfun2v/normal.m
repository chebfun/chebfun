function N = normal( F, varargin )
%NORMAL normal vector to a surface represented by a CHEBFUN2V.
%   N = NORMAL(F) returns a CHEBFUN2V representing the normal vector to the
%   surface F. The vector has the same magntiude as the surface's tangent vector
%
%   N = NORMAL(F,'unit') returns the unit normal vector, represented as a
%   CHEBFUN2V, to the surface F.
%
% See also CHEBFUN/NORMAL. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

N = cross( diff( F, 1, 2 ), diff( F, 1, 1 ) );

if ( nargin > 1 )
    if ( strcmpi(varargin{1}, 'unit') )
        if ( N.nComponents == 2 )
            r = roots( N );
            if ( ( ~isempty(r) ) || ( norm(N) == 0 ) )
                error('CHEBFUN:CHEBFUN2V:normal:zeroNormal1', ...
                    'Normal vector is zero.');
            end
            N = N ./ abs( N );
        elseif ( norm(N) == 0 )
            error('CHEBFUN:CHEBFUN2V:normal:zeroNormal2', ...
                'Normal vector is zero');
        end          
        
    else
        error('CHEBFUN:CHEBFUN2V:normal:badInput', ...
            'Second argument is not recognised.');
    end
end

end
