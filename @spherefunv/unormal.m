function N = unormal( dom )
%UNORMAL unit normal vector to the unit sphere.
%   N = UNORMAL( DOM ) returns a SPHEREFUNV representing the unit normal
%   vector to the unit sphere using the spherical coordinate domain
%   specified, i.e. DOM = [-pi -pi 0 pi], for co-latitude.
%
%   N = UNORMAL returns a SPHEREFUNV representing the unit normal
%   vector to the unit sphere, using the default spherical coordinates
%   domain, which is co-latitude.
%
%   See also TANGENTIAL, DOT, CURL

if ~exist( 'dom', 'var' ) || isempty( dom )
    dom = [-pi pi 0 pi];
end

if ( numel( dom ) ~= 4 )
    error('SPHEREFUN:SPHEREFUNV:normal:domain',...
        ['Domain for the unit vector is incorrect.  It should be a '...
        'double vector with 4 components.']);
end

% Normal vector to the unit sphere.
x = spherefun(@(x,y,z) x, dom);
y = spherefun(@(x,y,z) y, dom);
z = spherefun(@(x,y,z) z, dom);
N = spherefunv(x,y,z);

end
