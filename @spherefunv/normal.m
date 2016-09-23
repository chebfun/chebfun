function N = normal(F)
%NORMAL   Projection onto the normal vector to the unit sphere.
%   N = NORMAL(F) returns a SPHEREFUNV representing the projection of F
%   onto the normal vector to the sphere.
%
% See also SPHEREFUNV/TANGENT, SPHEREFUNV/DOT, SPHEREFUNV/CURL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If F is not a SPHEREFUNV then:
if ( ~isa(F, 'spherefunv') )
    error('SPHEREFUN:SPHEREFUNV:normal:type', 'Input must be a spherefunv.');
end

% If F is an empty SPHEREFUNV then:
if ( isempty(F) )
    N = spherefunv;
    return;
end

Fc = F.components;
dom = Fc{1}.domain;

% Normal vector to the unit sphere.
x = spherefun(@(x,y,z) x, dom);
y = spherefun(@(x,y,z) y, dom);
z = spherefun(@(x,y,z) z, dom);
N = spherefunv(x,y,z);

% If the spherefunv is empty just return the normal vector to the sphere.
if isempty(F)
    return;
end

% Dot F with the unit formal.
f = dot(F, N);

% Multiply by the normal.
N = N.*f;

end
