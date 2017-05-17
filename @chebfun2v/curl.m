function F = curl(F)
%CURL  curl of a CHEBFUN2V
%   S = CURL(F) returns the CHEBFUN2 of the curl of F. If F is a CHEBFUN2V with
%   two components then it returns the CHEBFUN2 representing
%         CURL(F) = F(2)_x - F(1)_y,
%   where F = (F(1),F(2)).  If F is a CHEBFUN2V with three components then it
%   returns the CHEBFUN2V representing the 3D curl operation.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.  

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
