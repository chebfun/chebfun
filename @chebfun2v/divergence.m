function f = divergence( F )
%DIVERGENCE the divergence of a chebfun2v.
%
% f = DIVERGENCE(F) returns the divergence of the chebfun2v i.e. 
% 
%  divergence(F) = F_x + F_y

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

% Note that divergence of a 3-vector is the same, because the functions are
% of two variables.

% Empty check: 
if ( isempty( F ) )
    f = chebfun2;
    return
end

Fc = F.components; 
f = diff(Fc(1), 1, 2) + diff( Fc(2), 1, 1 );  % divergence = f_x + f_y

end