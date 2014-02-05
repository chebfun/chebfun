function L = laplacian( F )
%LAPLACIAN Vector Laplacian of a chebfun2v.
%   LAPLACIAN(F) returns a chebfun2v representing the vector Laplacian of F.
%
% See also CHEBFUN2V/LAP.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty( F ) )
    L = chebfun2v;
    return
end

% laplacian = f_xx + f_yy
L = diff(F, 2, 1) + diff(F, 2, 2);   

end