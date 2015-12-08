function pass = testspherefun( )
% Main testing file, for now.  
pass(1) = all(test_constructor( )); 
pass(2) = all(test_feval( )); 
pass(3) = all(test_sum2( )); 
pass(4) = all(test_plus( ));
pass(5) = all(test_times( )); 
pass(6) = all(test_power( ));
pass(7) = all(test_abs( ));
pass(8) = all(test_diff( ));
pass(9) = all(test_laplacian( ));
end 




















function [x, y] = getPoints( m, n )

% x = linspace(-pi, pi, 2*n+1)';  x( end ) = [ ];
% GBW: You can't just remove the pole and keep everything else equally
% spaced between -pi/2 and 3*pi/2.  The issue is that you can't keep
% both point -pi/ and 3*pi/2.
% y = linspace(-pi/2, 3*pi/2, 2*m+1)'; y( m+1 ) = [ ];

x = trigpts(2*n,[-pi pi]);
% GBW: I believe we have to sample at equally spaced points shifted by h/2
% to not sample the poles and keep an even total number of points.
% y = trigpts(2*m,[-pi/2 3*pi/2]);
% y = trigpts(m,[0 pi]);
y = linspace(0,pi,m).';

% y = y+0.5*pi/m; % Shift y by h/2 to avoid the poles
end

function f = redefine_function_handle( f )
% nargin( f ) = 2, then we are already on the sphere, if nargin( f ) = 3,
% then do change of variables:

if ( nargin( f ) == 3 )
    % Wrap f so it can be evaluated in spherical coordinates
    f = @(lam, th) spherefun.sphf2cartf(f,lam,th,0);
%     % Double g up.
%     f = @(lam, th) sph2torus(f,lam,th);
end

end


% function fdf = sph2torus(f,lam,th)
% 
% fdf = real(f(lam,th));
% 
% id = th-pi/2 > 100*eps;
% 
% if ~isempty(id) && any(id(:))
%     fdf(id) = f(lam(id)-pi,pi-th(id));
% end
% 
% end
% 
% function fdf = sphf2cartf(f,lam,th)
% 
% x = cos(lam).*cos(th);
% y = sin(lam).*cos(th);
% z = sin(th);
% 
% fdf = f(x,y,z);
% 
% end


