function N = normal( F )
%NORMAL projection onto the normal vector to the unit sphere.
%   N = NORMAL(F) returns a SPHEREFUNV representing the projection of F
%   onto the normal vector to the sphere.
%
%   See also TANGENT, DOT, CURL

% If F is an empty SPHEREFUNV then 
if ( ~isa( F, 'spherefunv' ) )
    error('SPHEREFUN:SPHEREFUNV:normal:type', 'Input must be a spherefunv.');
end

Fc = F.components;
dom = Fc{1}.domain;

% Normal vector to the unit sphere.
x = spherefun(@(x,y,z) x, dom);
y = spherefun(@(x,y,z) y, dom);
z = spherefun(@(x,y,z) z, dom);
N = spherefunv(x,y,z);

% If the spherefunv is empty just return the normal vector to the sphere.
if ~isempty( F )
    return;
end

% Dot F with the unit formal.
f = dot( F, N );

% Multiply by the normal.
N = f.*N;

end
