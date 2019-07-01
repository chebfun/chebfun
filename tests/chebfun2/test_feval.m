function pass = test_feval(pref)
% Test feval

if ( nargin == 0)
    pref = chebfunpref;
end

tol = 100*pref.cheb2Prefs.chebfun2eps;

f = chebfun2(@(x,y) x);
pass(1) = (abs(f(pi/6,pi/12)-pi/6) < tol);

f = chebfun2(@(x,y) x, [-1 2 -pi/2 pi]);
pass(2) = (abs(f(0,0)) < 1e-14);
pass(3) = (abs(f(pi/6,pi/12)-pi/6) < tol);


f = chebfun2(@(x,y) y, [-1 2 -pi/2 pi]);
pass(4) = (abs(f(0,0)) < 1e-14);
pass(5) = (abs(f(pi/6,pi/12)-pi/12) < tol);


% some harder tests.
f = @(x,y) cos(x) + sin(x.*y);
g = chebfun2(f);

r = 0.126986816293506; s = 0.632359246225410; % two fixed random number in domain.
pass(6) = (abs(f(r,s) - g(r,s))<tol);

% Are we evaluating on arrays correctly
seedRNG(0)
r = rand(10,1);
s = rand(10,1);
[rr, ss]=meshgrid(r,s);
pass(7) = (norm((f(r,s) - g(r,s))) < tol);
pass(8) = (norm((f(rr,ss) - g(rr,ss))) < tol);  % on arrays as well.

% Does this work off [-1,1]^2
g = chebfun2(f,[-pi/6 pi/2 -pi/12 sqrt(3)]); % strange domain.
r = 0.126986816293506; s = 0.632359246225410; % two fixed random number in domain.
pass(9) = (abs(f(r,s) - g(r,s))<tol);

% Are we evaluating on arrays correctly
r = rand(10,1); s = rand(10,1); [rr, ss]=meshgrid(r,s);
pass(10) = (norm((f(r,s) - g(r,s)))<tol);
pass(11) = (norm((f(rr,ss) - g(rr,ss)))<2*tol); % on arrays as well.

% Evaluation at complex arguments.
f = chebfun2(@(x,y) x+1i*y);
pass(12) = ( norm( f(-1,-1) - (-1-1i) )  < tol );
pass(13) = ( norm( f(r,s) - (r + 1i*s) )  < tol );
pass(14) = (norm((f(rr,ss) - (rr + 1i*ss)))<tol); % on arrays as well.

% Evaluation at complex arguments, different syntax.
f = chebfun2(@(z) z);
pass(15) = ( norm( f(-1-1i) - (-1-1i) )  < tol );
pass(16) = ( norm( f(r+1i*s) - (r + 1i*s) )  < tol );
pass(17) = (norm((feval(f,rr,ss) - (rr + 1i*ss)))<tol); % on arrays as well.

% Evaluation at transposed meshgrid
op = @(x,y) cos(pi*(x+y))+y.*sin(pi*x);
[yy,xx] = meshgrid(linspace(-1,1,1001));
f = chebfun2( op );
pass(18) = (norm( feval(f,xx,yy) - op(xx,yy) ) < 100*tol);

% Check sizes:
n = 10;
f = chebfun2(@(x,y) cos(x.*y));
x = ones(1,n);
[xx, yy] = meshgrid(x);
pass(19) = ( all( size( feval(f, 1, 1) ) == [1 1] ) );
pass(20) = ( all( size( feval(f, [1 1], [1 1]) ) == [1 2] ) );
pass(21) = ( all( size( feval(f, [1;1], [1;1]) ) == [2 1] ) );
pass(22) = ( all( size( feval(f, [1 1;1 1], [1 1; 1 1]) ) == [2 2] ) );
pass(23) = ( all( size( feval(f, x, x) ) == size(x) ) );
pass(24) = ( all( size( feval(f, xx, yy)  ) == [n n] ) );

end
