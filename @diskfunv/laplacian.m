function L = laplacian( F )
%LAPLACIAN Vector Laplacian of a DISKFUNV.
%   LAPLACIAN(F) returns a DISKFUNV representing the vector Laplacian of F.
% 
% See also DISKFUN/LAPLACIAN

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) )
    L = diskfunv();
    return
end

% The vector LAPLACIAN of a DISKFUNV is equal to F_xx + F_yy: 
L = diff(F, 1, 2) + diff(F, 2, 2);   

end
