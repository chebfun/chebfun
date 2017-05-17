function vol = domainvolume(f)
%DOMAINVOLUME   Volume of the domain of a CHEBFUN3 object.
%   DOMAINVOLUME(F) returns the volume of the domain of F.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    vol = 0;
else
    dom = f.domain; 
    vol = diff(dom(1:2)) * diff(dom(3:4)) * diff(dom(5:6));
end

end