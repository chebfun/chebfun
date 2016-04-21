function pass = test_guide( pref )
% Test Chebfun2 guide commands. This is not exclusive by checks the main 
% commands. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb2Prefs.chebfun2eps; 

% commands from guide1: 
f = chebfun2(@(x,y) cos(x.*y));

pass(1) = abs( sum2(f) - 3.784332281468732 ) < tol; 

x = chebfun2(@(x,y) x, [-2 3 -4 4]); 
y = chebfun2(@(x,y) y, [-2 3 -4 4]);   

op = @(x,y) 1./( 2 + cos(.25 + x.^2.*y + y.^2) );
g = ( 2 + cos(.25 + x.^2.*y + y.^2) );
f = 1./g;

x = linspace(-2,3); 
y = linspace(-4,4); 
[xx,yy] = meshgrid(x,y); 
pass(2) = ( norm( f(xx,yy) - op(xx,yy), inf) ) < 200*tol; 

f = chebfun2(@(z) sin(z));   
pass(3) = ( abs( f(1+1i) -  sin(1+1i) ) ) < tol;

f = chebfun2(@(z) sin(z)-sinh(z),2*pi*[-1 1 -1 1]);
x = linspace(-2*pi,2*pi);  
[xx,yy] = meshgrid(x); 
pass(4) = ( norm( f(xx+1i*yy) - f(xx,yy), inf)) < tol;
pass(5) = ( norm( f(xx+1i*yy) - sin(xx+1i*yy)+sinh(xx+1i*yy), inf)) < 1e3*tol;

% commands from guide2: 
f = chebfun2(@(x,y) sin(10*x.*y), [0 pi/4 0 3]);
pass(6) = ( norm( sum(f) - chebfun(@(x) sin(15*x).^2./x/5,[0 pi/4])' )) < tol;
pass(7) = ( norm( sum(f,2) - chebfun(@(y) sin(5*pi*y/4).^2./y/5,[0 3]) )) < tol;
pass(8) = ( abs( sum2(f) - sum(sum(f)) )) < tol;

F = @(x,y) exp(-(x.^2 + y.^2 + cos(4*x.*y))); 
I1 = quad2d(F,-1,1,-1,1,'AbsTol',tol);
I2 = sum(sum(chebfun2(F)));
I3 = sum2(chebfun2(F)); 
pass(9) = ( abs(I1 - I2)) < tol;
pass(10) = ( abs(I2 - I3)) < tol;
pass(11) = ( abs(I1 - I3)) < tol;

f = chebfun2( 'exp(-(x.^2 + y.^2 +4*x.*y))' );
pass(12) = ( abs( norm(f) - sqrt(sum2(f.^2)) ) < tol );

f = chebfun2( @(x,y) exp(-1./( sin(x.*y) + x ).^2) );
g = chebfun2( @(x,y) cos(exp(-1./( sin(x.*y) + x ).^2)) );
h = chebfun2( @(x,y) (exp(-1./( sin(x.*y) + x ).^2)).^5 );
pass(13) = ( norm( cos(f) - g ) ) < tol;
pass(14) = ( norm( f.^5 - h )) < tol;

runge = chebfun2( @(x,y) 1./( .01 + x.^2 + y.^2 )) ;
pass(15) = ( abs(mean2(runge) - 3.796119578934828)) < tol;
pass(16) = ( abs(mean(mean(runge)) - 3.796119578934828)) < tol;

f = chebfun2(@(x,y) exp(-(x.^2 + 3*x.*y+y.^2) ));
pass(17) = ( norm( cumsum(cumsum(f),2) - cumsum2(f) )) < tol;

f = chebfun2(@(x,y) cos(10*x.*y.^2) + exp(-x.^2)); % chebfun2 object
C = chebfun(@(t) t.*exp(10i*t),[0 1]);             % spiral curve

pass(18) = ( abs(sum(f(C)) - 1.613596461872283)) < tol;

% commands from guide4: 
f = chebfun2(@(x,y) sin(x+1i*y));   % a holomorphic function
u = real(f); v = imag(f);           % real and imaginary parts

pass(19) = ( norm(diff(u) - (-diff(v,1,2)))) < tol;
pass(20) = ( norm(diff(u,1,2) - diff(v)) ) < tol;        % Do the Cauchy-Riemann eqns hold?

d = [0 1 0 2];
F = chebfun2v(@(x,y) sin(x.*y), @(x,y) cos(y),d);  % calling the constructor
f = chebfun2(@(x,y) sin(x.*y),d); g = chebfun2(@(x,y) cos(y),d);
G = [f;g];
plaw = abs((2*norm(F)^2 + 2*norm(G)^2) - (norm(F+G)^2 + norm(F-G)^2));
pass(21) = ( plaw < tol ); 

f = chebfun2(@(x,y) cos(10*x.*y.^2) + exp(-x.^2));  % chebfun2
F = gradient(f);                                    % gradient (chebfun2v)
C = chebfun(@(t) t.*exp(10i*t),[0 1]);              % spiral curve
v = integral(F,C);ends = f(cos(10),sin(10))-f(0,0); % line integral
pass(22) = ( abs(v-ends) ) < tol;

% commands from guide6: 
r1 = 1; r2 = 1/3;   % inner and outer radius
d = [0 2*pi 0 2*pi];
u = chebfun2(@(u,v) u,d);
v = chebfun2(@(u,v) v,d);
F = [-(r1+r2*cos(v)).*sin(u); (r1+r2*cos(v)).*cos(u); r2*sin(v)];  % torus
G = F./3;  % full 3D divergence of G is 1 because F = [x;y;z]. 
pass(23) = ( abs( integral2(dot(G,normal(F))) - 2*pi^2*r1*r2.^2 ) < tol );
