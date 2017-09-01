function T = tangent( F )
%TANGENT   Projection onto the tangent plane of the unit sphere.
%   T = TANGENT(F) returns a SPHEREFUNV representing the projection of F
%   onto the tangent plane to the unit sphere. The result is a vector
%   field in the tangent space of the unit sphere.
%
%   See also SPHEREFUNV/NORMAL, SPHEREFUNV/DOT, SPHEREFUNV/CURL

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% If the spherefunv is empty just return an empty spherefunv.
if isempty(F)
    T = F;
    return;
end

% The projection is just F - normal(F).
T = F - normal(F);

end

