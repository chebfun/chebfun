function T = tangent( F )
%TANGENT   projection onto the tangent plane of the unit sphere.
%   T = TANGENT(F) returns a SPHEREFUNV representing the projection of F
%   onto the tangent plane to the unit sphere.  The result is a vector
%   field in the tangent space of the unit sphere.
%
%   See also NORMAL, DOT, CURL

% If the spherefunv is empty just return an empty spherefunv.
if isempty( F )
    T = F;
    return;
end

% The projection is just F - normal(F).
T = F - normal( F );

end

