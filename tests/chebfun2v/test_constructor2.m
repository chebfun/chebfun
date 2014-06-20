function pass = test_constructor2( pref )
% Test the Chebfun2v constructor when performing simple arithmetic
% operations.

if ( nargin < 1 )
    pref = chebfunpref;
end

tol = 100*pref.eps;

% check the vectorize flag: 
f1 = @(x,y) x.*y; 
f2 = @(x,y) x*y; 
H1 = chebfun2v(f1, f1);
H2 = chebfun2v(f2, f2, 'vectorize');
pass(1) = norm( H1 - H2 ) < tol; 

f1 = @(x,y) x.*y; 
f2 = @(x,y) x*y; 
H1 = chebfun2v(f1, f1, [-2 3 -1 0]);
H2 = chebfun2v(f2, f2, 'vectorize', [-2 3 -1 0]);
pass(2) = norm( H1 - H2 ) < tol; 

f1 = @(x,y) x.*y; 
f2 = @(x,y) x*y; 
H1 = chebfun2v(f1, f1, f1);
H2 = chebfun2v(f2, f2, f2, 'vectorize');
pass(3) = norm( H1 - H2 ) < tol; 

f1 = @(x,y) x.*y; 
f2 = @(x,y) x*y; 
H1 = chebfun2v(f1, f1, f1, [-2 3 -1 0]);
H2 = chebfun2v(f2, f2, f2, 'vectorize', [-2 3 -1 0]);
pass(4) = norm( H1 - H2 ) < tol; 

% Test the constructor with a surface example: 
u = chebfun2(@(u,v) u, [0 2*pi -1 1]);
v = chebfun2(@(u,v) v, [0 2*pi -1 1]);
x = (1+0.5*v.*cos(u/2)).*cos(u);
y = (1+0.5*v.*cos(u/2)).*sin(u);
z = 0.5*v.*sin(u/2);
r = [x;y;z];
ru = diff(r,1,1);
rv = diff(r,1,2);
pass(5) = ( norm(ru'*rv,inf) < 10*tol ); 

V = [sin(5*u);cos(5*v);0];    
pass(6) = (V.nComponents == 3); 

end