function f = divergence( F )
%DIVERGENCE   Divergence of a ADCHEBFUN2.
%   DIVERGENCE(F) returns the divergence of the ADCHEBFUN2V i.e.,
%   divergence(F) = F{1}_x + F{2}_y.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information. 

% Note that divergence of a 3-vector is the same, because the functions are
% of two variables.

% divergence = f_x + f_y
f = diff( F(1), 1, 2 ) + diff( F(2), 1, 1 );  

end