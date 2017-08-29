function [u, v] = helmholtzdecomp( f ) 
%HELMHOLTZDECOMP   Return the Helmholtz decomposition of spherefunv.
%
% [U, V] = HELMHOLTZDECOMP( F ) computes the Helmholtz decomposition
% of the spherefunv F, i.e., 
%        F   =     GRAD( U )     +    CURL( V ) 
% where U and V are spherefun objects. F needs to be a vector field that is
% tangential to the surface of the sphere. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    u = spherefun; 
    v = spherefun; 
    return
end

tangentf = tangent( f );
tol1 = 100*vscale(f.components{1})*chebfunpref().cheb2Prefs.chebfun2eps;
tol2 = 100*vscale(f.components{2})*chebfunpref().cheb2Prefs.chebfun2eps;
tol3 = 100*vscale(f.components{3})*chebfunpref().cheb2Prefs.chebfun2eps;
if ( norm( f - tangentf ) ) > tol1+tol2+tol3 
    warning( 'SPHEREFUNV:HELMHOLTZDECOMPOSITON:TANGENT', ...
        ['The vector field needs to be tangent to the surface of the',...
        'sphere, taking f = tangent(f).']) 
end 

% Divergence of curl is zero => div(f) = div(grad(u)).  Also, 
% div(grad(u)) = lap(u): 
divf = div(tangentf); 
[m, n] = length( divf ); 
% Just pick a sufficiently large discretization size. 
u = spherefun.poisson( divf, 0, max(2*m,50), max(2*n,50) ); 

% Vorticity of grad is zero on the surface of the spherefun => 
% vort( f ) = vort( curl( v ) ). Also, vort( curl( v ) ) = lap( v ) on the
% surface of the sphere: 
vortf = vort(tangentf); 
[m,n] = length( vortf ); 
v = spherefun.poisson( vortf, 0, max(2*m,50), max(2*n,50) );

end