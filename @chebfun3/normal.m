function N = normal(F, varargin)
%NORMAL normal vector to an implicit surface represented by a CHEBFUN3.
%
%   N = NORMAL(F) returns a CHEBFUN3V representing the normal vector to the
%   surface F. The vector has the same magntiude as the surface's tangent 
%   vector.
%
%   N = NORMAL(F, 'unit') returns the unit normal vector, represented as a
%   CHEBFUN3V, to the surface F.
%
% See also CHEBFUN/NORMAL and CHEBFUN2V/NORMAL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

N = grad(F);

if ( nargin > 1 )
    if ( strcmpi(varargin{1}, 'unit') )
        nrm = norm(N);
        if ( nrm == 0 )
            error('CHEBFUN:CHEBFUN2V:normal:zeroNormal2', ...
                'Normal vector is zero');
        else
            N = N ./ nrm;
        end
    end
end

end