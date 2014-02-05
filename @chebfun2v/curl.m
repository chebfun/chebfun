function F = curl(F)
%CURL  curl of a chebfun2v
%   S = CURL(F) returns the chebfun2 of the curl of F. If F is a chebfun2v with
%   two components then it returns the chebfun2 representing
%         CURL(F) = F(2)_x - F(1)_y,
%   where F = (F(1),F(2)).  If F is a chebfun2v with three components then it
%   returns the chebfun2v representing the 3D curl operation.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.  

% Empty check: 
if ( isempty( F ) )
    F = chebfun2v;
    return
end

Fc = F.components; 

if ( F.nComponents == 2 )  % 2D curl 
    F = diff(Fc{2}, 1, 2) - diff(Fc{1}, 1, 1);
else                       % 3D curl
    F = [ diff(Fc{3},1,1) ; -diff(Fc{3},1,2) ;...
          diff(Fc{2},1,2) - diff(Fc{1},1,1) ];
end

end