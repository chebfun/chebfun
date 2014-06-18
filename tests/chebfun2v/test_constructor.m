function pass = test_constructor( pref ) 
% Test the Chebfun2v constructor when performing simple arithmetic
% operations. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = 100*pref.eps; 

D = [-1 1 -1 1; -2 3 0 1];

for kk = 1:size(D,1)
    d = D(kk,:);
% Check the constructor works with lots of different syntax: 
f1 = @(x,y) cos(x.*y); 
f2 = @(x,y) sin(x) + cos(y); 
f3 = @(x,y) exp(x) - exp(-y) + x; 
g1 = chebfun2(f1, d);
g2 = chebfun2(f2, d); 
g3 = chebfun2(f3, d);

H1 = chebfun2v(f1, f2, d);
H2 = chebfun2v(f1, g2, d);
H3 = chebfun2v(g1, f2, d);
H4 = chebfun2v(g1, g2, d);

pass(1) = norm( H1 - H2 ) < tol; 
pass(2) = norm( H1 - H3 ) < tol; 
pass(3) = norm( H1 - H4 ) < tol; 

% Now try the same thing for three components: 
H1 = chebfun2v(f1, f2, f3, d);
H2 = chebfun2v(f1, f2, g3, d);
H3 = chebfun2v(f1, g2, f3, d);
H4 = chebfun2v(g1, f2, f3, d);
H5 = chebfun2v(g1, g2, f3, d);
H6 = chebfun2v(g1, f2, g3, d);
H7 = chebfun2v(f1, g2, g3, d);
H8 = chebfun2v(g1, g2, g3, d);

pass(4) = norm( H1 - H2 ) < tol; 
pass(5) = norm( H1 - H3 ) < tol; 
pass(6) = norm( H1 - H4 ) < tol; 
pass(7) = norm( H1 - H5 ) < tol; 
pass(8) = norm( H1 - H6 ) < tol; 
pass(9) = norm( H1 - H7 ) < tol; 
pass(10) = norm( H1 - H8 ) < tol; 

end

% check the vectorize flag: 
f1 = @(x,y) x.*y; 
f2 = @(x,y) x*y; 
H1 = chebfun2v(f1, f1);
H2 = chebfun2v(f2, f2, 'vectorize');
pass(11) = norm( H1 - H2 ) < tol; 

f1 = @(x,y) x.*y; 
f2 = @(x,y) x*y; 
H1 = chebfun2v(f1, f1, [-2 3 -1 0]);
H2 = chebfun2v(f2, f2, 'vectorize', [-2 3 -1 0]);
pass(12) = norm( H1 - H2 ) < tol; 

f1 = @(x,y) x.*y; 
f2 = @(x,y) x*y; 
H1 = chebfun2v(f1, f1, f1);
H2 = chebfun2v(f2, f2, f2, 'vectorize');
pass(13) = norm( H1 - H2 ) < tol; 

f1 = @(x,y) x.*y; 
f2 = @(x,y) x*y; 
H1 = chebfun2v(f1, f1, f1, [-2 3 -1 0]);
H2 = chebfun2v(f2, f2, f2, 'vectorize', [-2 3 -1 0]);
pass(14) = norm( H1 - H2 ) < tol; 

% Test the constructor with a surface example: 
u = chebfun2(@(u,v) u, [0 2*pi -1 1]);
v = chebfun2(@(u,v) v, [0 2*pi -1 1]);
x = (1+0.5*v.*cos(u/2)).*cos(u);
y = (1+0.5*v.*cos(u/2)).*sin(u);
z = 0.5*v.*sin(u/2);
r = [x;y;z];
n = normal(r);
ru = diff(r,1,1);
rv = diff(r,1,2);
pass(15) = ( norm(ru'*rv,inf) < tol ); 

R1 = ru/sqrt(ru'*ru);
R2 = rv/sqrt(rv'*rv);
V = [sin(5*u);cos(5*v);0];    
% Can we construct PV?:
pass(16) = 1; 
try 
    PV = (R1'*V)*R1 + (R2'*V)*R2;
catch
    pass(16) = 0; 
end

end