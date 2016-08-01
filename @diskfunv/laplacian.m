function L = laplacian( F )
%LAPLACIAN Vector Laplacian of a DISKFUNV.
%   LAPLACIAN(F) returns a DISKFUNV representing the vector Laplacian of F.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) )
    L = diskfunv;
    return
end

% laplacian = f_xx + f_yy
L = diff(F, 1, 2) + diff(F, 2, 2);   

end
