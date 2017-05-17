function L = laplacian( F )
%LAPLACIAN Vector Laplacian of a CHEBFUN2V.
%   LAPLACIAN(F) returns a CHEBFUN2V representing the vector Laplacian of F.
%
% See also CHEBFUN2V/LAP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) )
    L = chebfun2v;
    return
end

% laplacian = f_xx + f_yy
L = diff(F, 2, 1) + diff(F, 2, 2);   

end
