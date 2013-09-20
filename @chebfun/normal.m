function n = normal(c,unit)
%NORMAL   Normal to a complex-valued CHEBFUN.
%   N = NORMAL(C) returns the normal vector to the curve C as a quasi-matrix
%   with two columns. The vector has the same magntiude as the curve's tangent
%   vector.
%
%   N = NORMAL(C, 'unit') returns the unit normal vector to the curve C. N is a
%   quasi-matrix with two columns.
% 
% See also CHEBFUN2V/NORMAL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TOD0]: This method requires HORZCAT() and quasimatrices.

dc = diff(c);
n = -1i*dc; 

if ( nargin > 1 ) 
    if ( strcmpi(unit, 'unit') )
        n = n./norm(n);
    else
        error('CHEBFUN:NORMAL', 'Second argument is not recognised.');
    end
end

% Return a quasi-matrix. 
n = [real(n), imag(n)];  
 
end