function pass = test_roots( )
% Test the chebfun2/roots command.

tol = 1e-12;  % rank 1 curves
tol2 = 1e-8;  % rank >1 curves, disjoint
tol3 = 1e-4;  % rank >1 curves, noncircular
tol4 = 1e-2;  % intersecting curves, not properly disentangled

% First we check separable functions, i.e., rank 1:

f = chebfun2( @(x,y) 1/2-y);
r = roots( f );
exact = chebfun( @(x) x + 1i/2 );
pass(1) = ( norm( r - exact ) < tol );

f = chebfun2( @(x,y) 1/2-x);
r = roots( f );
exact = chebfun( @(x) 1i*x + 1/2 );
pass(2) = ( norm( r - exact ) < tol );

f = chebfun2( @(x,y) x.*y);
r = roots( f );
exact = chebfun( @(x) [x, 1i*x]  );
pass(3) = ( norm( r - exact ) < tol );

f = chebfun2( @(x,y) cos(5*pi*x));
s = roots( squeeze( f ) );
r = roots( f );
exact = chebfun;
for j = 1:numel(s)
    exact = [ exact chebfun(@(x) s(j) + 1i*x )];
end
pass(4) = ( norm( r - exact ) < tol );

f = chebfun2( @(x,y) x.*(y-1/2), [-2 2 -3 5]);
r = roots( f );
exact = [chebfun( @(x) 2*x + 1i/2 ) chebfun(@(x) 1i*(4*(x+1)-3) )];
pass(5) = ( norm( r - exact ) < tol );

% Next we move to functions of rank > 1

f = chebfun2( @(x,y) x.^2 + y.^2 - 1/4 );
c = roots(f);
arclength = sum(abs(diff(c)));
err = arclength - pi;
pass(6) = abs(err) < tol2;
%disp([int2str(pass(6)) '  circle  ' int2str(length(c))])

area = abs(sum(real(c).*diff(imag(c))));
err = area - pi/4;
pass(7) = abs(err) < tol2;
%disp([int2str(pass(7)) '  area of circle  ' int2str(length(c))])

f = chebfun2( @(x,y) x.^2 + (10*y).^2 - 1/4 );
c = roots(f);
arclength = sum(abs(diff(c)));
err = arclength - 2.031987090050447;
pass(8) = abs(err) < tol3;
%disp([int2str(pass(8)) '  ellipse  ' int2str(length(c))])

area = abs(sum(real(c).*diff(imag(c))));
err = area - pi/40;
pass(9) = abs(err) < tol3;
%disp([int2str(pass(9)) '  area of ellipse  ' int2str(length(c))])

f = chebfun2( @(x,y) (x-1).^2 + y.^2 - 1/4 );
c = roots(f);
arclength = sum(abs(diff(c)));
err = arclength - pi/2;
pass(10) = abs(err) < tol2;
%disp([int2str(pass(10)) '  semicircle at boundary  ' int2str(length(c))])

f = chebfun2( @(x,y) (x-1).^2 + (y+1).^2 - 1/4 );
c = roots(f);
arclength = sum(abs(diff(c)));
err = arclength - pi/4;
pass(11) = abs(err) < tol2;
%disp([int2str(pass(11)) '  quarter-circle at corner  ' int2str(length(c))])

r = 1/pi;  % this radius gives circumference 2
x = chebfun2(@(x,y) x);
y = chebfun2(@(x,y) y);

%% One circle
f = x.^2 + y.^2 - r^2;      
c = roots(f); arclength = sum(abs(diff(c)));
err = arclength - 2;
pass(12) = abs(err) < tol2;
%disp([int2str(pass(12)) '  another circle  ' int2str(length(c))])

%% Two circles
f = x.^2 + y.^2 - r^2;      
g = (x-.6).^2 + (y-.5).^2 - r^2;      
c = roots(f.*g); arclength = sum(abs(diff(c)));
err = arclength - [2 2];
pass(13) = norm(err,inf) < tol2;
%disp([int2str(pass(13)) '  two circles  ' int2str(length(c))])

%% Intersecting circles
f = x.^2 + y.^2 - r^2;      
g = (x-.2).^2 + (y-.5).^2 - r^2;      
c = roots(f.*g); arclength = sum(abs(diff(c)));
err = sum(arclength) - 4;
pass(14) = abs(err) < tol4;
%disp([int2str(pass(14)) '  intersecting circles  ' int2str(length(c))])

%% A circle and a line
f = x.^2 + y.^2 - r^2;      
w = exp(x)-1.2;
c = roots(f.*w); arclength = sum(abs(diff(c)));
err = sum(arclength) - 4;
pass(15) = abs(err) < tol4;
%disp([int2str(pass(15)) '  intersecting circle and line  ' int2str(length(c))])

%% Concentric circles (well separated)
f = chebfun2(@(x,y) (x.^2 + y.^2 - .25).*(x.^2+y.^2-0.300));
c = roots(f); arclength = sum(abs(diff(c)));
err = arclength - pi*[sqrt(.300/.250) 1];
pass(16) = norm(err,inf) < tol2;
%disp([int2str(pass(16)) '  concentric circles  ' int2str(length(c))])

%% Concentric circles (less well separated)
f = chebfun2(@(x,y) (x.^2 + y.^2 - .25).*(x.^2+y.^2-0.260));
c = roots(f); arclength = sum(abs(diff(c)));
err = arclength - pi*[sqrt(.260/.250) 1];
pass(17) = norm(err,inf) < tol2;
%disp([int2str(pass(17)) '  close concentric circles  ' int2str(length(c))])

f = chebfun2( @(x,y) (x.^2+y.^2-1/4).*((x-1.1).^2+(y-3).^2-1/4),[-1.3 1.1 -.9 3]);
c = roots(f);
arclength = sum(abs(diff(c)));
err = sum(arclength) - 1.25*pi;
pass(18) = abs(err) < tol2;
%disp([int2str(pass(18)) '  different domain  ' int2str(length(c))])

f = chebfun2( @(x,y) 1e100*(x.^2 + y.^2 - 1/4) );
c = roots(f);
arclength = sum(abs(diff(c)));
err = arclength - pi;
pass(19) = abs(err) < tol2;
%disp([int2str(pass(19)) '  vertical scaling 1  ' int2str(length(c))])

f = chebfun2( @(x,y) 1e100*(x.^2 + y.^2 - 1/4) );
area = abs(sum(real(c).*diff(imag(c))));
err = area - pi/4;
pass(20) = abs(err) < tol2;
%disp([int2str(pass(20)) '  vertical scaling 2  ' int2str(length(c))])

f = chebfun2( @(x,y) y-sin(x),[-pi pi,-1.1 1.3]);
c = roots(f);
arclength = sum(abs(diff(c)));
err = arclength - 7.640395578055424035;
pass(21) = abs(err) < tol2;

f = chebfun2( @(x,y) x-tanh(y));
c = roots(f);
arclength = sum(abs(diff(c)));
err = arclength - 2.53182746833076;
pass(22) = abs(err) < tol2;

f = chebfun2( @(x,y) y-x.^2);
c = roots(f);
arclength = sum(abs(diff(c)));
err = arclength - 2.957885715089195;
pass(23) = abs(err) < tol2;

f = chebfun2( @(x,y) y.^2-x.^2-.1,[-1 1 -1.5 1.5]);  % hyperbola
c = roots(f);
arclength = sum(abs(diff(c)));
err = arclength - [1 1]*2.51829399795912;
pass(24) = norm(err,inf) < tol2;

f = chebfun2( @(x,y) y.^2-x.^2);  % cross
c = roots(f);
arclength = sum(abs(diff(c)));
err = arclength - [1 1]*sqrt(8);
pass(25) = norm(err,inf) < tol4;

end
