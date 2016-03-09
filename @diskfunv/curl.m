function F = curl(F)
%CURL  curl of a DISKFUNV
%   S = CURL(F) returns the DISKFUN of the curl of F. If F is a DISKFUNV with
%   two components then it returns the DISKFUN representing
%         CURL(F) = F(2)_x - F(1)_y,
%   where F = (F(1),F(2)).  If F is a DISKFUNV with three components then it
%   returns the DISKFUNV representing the 3D curl operation.


% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.  

% Empty check: 
if ( isempty( F ) )
    F = diskfunv;
    return
end

Fc = F.components; 

%if ( F.nComponents == 2 )  % 2D curl 
    F = diff(Fc{2}, 1, 1) - diff(Fc{1}, 2, 1);
%else                       % 3D curl
    %F = [ diff(Fc{3},2,1) ; -diff(Fc{3},1,1) ;...
         % diff(Fc{2},1,1) - diff(Fc{1},2,1) ];
 %   F = diskfunv( diff(Fc{3}, 2), -diff(Fc{3},1), diff(Fc{2},1) - diff(Fc{1},2));
end

end
