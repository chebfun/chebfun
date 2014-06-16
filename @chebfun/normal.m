function n = normal(c, unit)
%NORMAL   Normal to a complex-valued CHEBFUN.
%   N = NORMAL(C) returns the normal vector to the curve C as a CHEBFUN with two
%   columns. The vector has the same magntiude as the curve's tangent vector.
%
%   N = NORMAL(C, 'unit') returns the unit normal vector to the curve C. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]:  Reconsider how this function should behave in the presence of cusps
% once we have singfun in place.

n = -1i*diff(c); 

if ( nargin > 1 ) 
    if ( strcmpi(unit, 'unit') )
        nrmn = norm(n);
        if ( nrmn == 0 )
            error('CHEBFUN:CHEBFUN:normal:zero', 'Normal vector is zero.'); 
        else
            n = n./nrmn;
        end
    else
        error('CHEBFUN:CHEBFUN:normal:args', ...
            'Second argument is not recognised.');
    end
end

% Return an array-valued CHEBFUN:
n = [real(n), imag(n)];  
 
end
