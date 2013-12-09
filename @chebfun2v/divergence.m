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

components = F.components; 
f = diff(components(1),1,2) + diff(components(2),1,1);  % divergence.

end