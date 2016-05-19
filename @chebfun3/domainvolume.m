function area = domainvolume(f)
%DOMAINVOLUME   Volume of the domain of a CHEBFUN3 object
%   DOMAINVOLUME(F) returns the volume of the topological domain of F.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    area = 0;
else
    dom = f.domain; 
    area = diff(dom(1:2)) * diff(dom(3:4)) * diff(dom(5:6));
end

end