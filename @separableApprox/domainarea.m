function area = domainarea( f )
%DOMAINAREA    Area of the domain of f
%
%   DOMAINAREA(F) returns the area of the topological domain of f.
%

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) )
    area = 0;
else
    dom = f.domain; 
    area = diff( dom(1:2) )*diff( dom(3:4) );   
end

end
