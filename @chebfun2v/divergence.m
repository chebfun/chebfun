function f = divergence( F )
%DIVERGENCE   Divergence of a CHEBFUN2V.
%   DIVERGENCE(F) returns the divergence of the CHEBFUN2V i.e.,
%       divergence(F) = F_x + F_y.
%
% See also DIV.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Note that divergence of a 3-vector is the same, because the functions are
% of two variables.

% Empty check: 
if ( isempty( F ) )
    f = chebfun2;
    return
end

Fc = F.components; 
% divergence = f_x + f_y
f = diff( Fc{1}, 1, 2 ) + diff( Fc{2}, 1, 1 );  

end
