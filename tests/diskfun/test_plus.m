function pass = test_plus( ) 
% Test diskfun plus() command 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f1 = @(x,y) sin(pi*x.*y);  % Strictly even/pi-periodic
f2 = @(x,y) sin(pi*x);  % Strictly odd/anti-periodic
g1 = diskfun(f1);
g2 = diskfun(f2);
gplus = g1 + g2;
fplus = @(t, r) diskfun.pol2cartf(@(x,y) f1(x, y) + ...
    f2(x, y), t, r);
t = pi*(2*rand-1); 
r = rand; 
pass(1) = abs(feval(gplus, t, r, 'polar') - fplus(t, r)) < tol;

f1 = @(x,y) sin(pi*x.*y) + sin(pi*x);   % Mixed symmetric terms
f2 = @(x,y) sin(pi*(x-.35).*(y-2.3)) + sin(pi*(x-.35)) ;  % Mixed symmetric terms
g1 = diskfun(f1);
g2 = diskfun(f2);
gplus = g1 + g2;
fplus = @(t, r) diskfun.pol2cartf( @(x,y) f1(x, y) + f2(x, y), t, r);
t = pi*(2*rand-1); 
r = rand;
pass(2) = abs(feval(gplus, t, r, 'polar') - fplus(t, r)) < tol; 

% Check that compression is working: 
f = diskfun(@(x,y) x.^2 + y.^2);
r = rank(f);
g = f; 
for k = 1:10
    g = g + f;
end 
pass(3) = norm(g - 11*f, inf) < vscale(g)*tol;
pass(4) = (rank(g) - r ) == 0;

% Check what happens with cancellation errors: 
f = diskfun(@(x,y) sin(x.*y)); 
g = 2*f;
pass(5) = ( norm(g - f - f, inf) < tol );

% Check for the case where f and g are non-zero at the poles bug 
% f+g is zero at the poles.
f = diskfun(@(x,y) x+1);
g = diskfun(@(x,y) cos(y));
h = f - g;
pass(6) = h.nonZeroPoles == 0;

% Check for the case where f non-zero at the poles but g is zero at the
% poles.
f = diskfun(@(x,y) cos(x));
g = diskfun(@(x,y) sin(y));
h = f + g;
pass(7) = h.nonZeroPoles == 1;
h = g + f;
pass(8) = h.nonZeroPoles == 1;

end