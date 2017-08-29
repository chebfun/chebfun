function F = curl(F)
%CURL  curl of a DISKFUNV
%   S = CURL(F) returns the DISKFUN of the curl of F. If F is a DISKFUNV 
%   then it returns the DISKFUN representing
%        S = CURL(F) = F(2)_x - F(1)_y,
%   where F = (F(1),F(2)).  

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.  

% Empty check: 
if ( isempty( F ) )
    F = diskfunv;
    return
end

% Extract the components of a DISKFUNV:
Fc = F.components; 

% Formula for 2D curl: 
F = diff(Fc{2}, 1, 1) - diff(Fc{1}, 2, 1);

end