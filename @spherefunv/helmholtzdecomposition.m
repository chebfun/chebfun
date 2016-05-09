function [u, v] = helmholtzdecomposition( f ) 
%HELMHOLTZDECOMPOSITION   Return the Helmholtz decomposition of spherefunv.
%
% [U, V] = HELMHOLTZDECOMPOSITION( F ) computes the Helmholtz decomposition
% of the spherefunv F, i.e., 
%        F   =     GRAD( U )     +    CURL( V ) 
% where U and V are spherefun objects. F needs to be a vector field that is
% tangential to the surface of the sphere. 
% 
% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

tangentf = tangent( f );

if ( norm( f - tangentf ) ) > tol 
    warning( 'SPHEREFUNV:HELMHOLTZDECOMPOSITON:TANGENT', ...
        ['The vector field needs to be tangent to the surface of the',...
        'sphere, taking f = tangent(f).']) 
end 

% Divergence of curl is zero => div(f) = div(grad(u)).  Also, 
% div(grad(u)) = lap(u): 
divf = div(f); 
[m,n] = length( divf ); 
u = spherefun.poisson( divf, 0, m, n ); 

% Vorticity of grad is zero on the surface of the spherefun => 
% vort( f ) = vort( curl( v ) ). Also, vort( curl( v ) ) = lap( v ) on the
% surface of the sphere: 
vortf = vort(f); 
[m,n] = length( vortf ); 
v = spherefun.poisson( vortf, 0, m, n );

end