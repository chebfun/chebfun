function pass = test_roots( pref ) 
% Check that roots works for a diskfun.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e5*pref.cheb2Prefs.chebfun2eps;

% More extensive testing is needed.  Here we just test that errors aren't
% thrown and the form of the output is reasonable.

f = diskfun(@(x,y) x+y );
r = roots(f);
if ~isempty(r)
    pass(1) = 1;
else
    pass(1) = 0;
end

% No zero contour exists
f = diskfun(@(x,y) 2+x.^2);
r = roots(f);
if isempty(r)
    pass(2) = 1;
else
    pass(2) = 0;
end

% Zero contour at  x.^2+y.^2=.5^2
f = diskfun(@(t,r) r.^2-.5^2, 'polar');
r = roots(f);

x = chebpts(257); 
pass(3) = (norm(.5*exp(1i*x*pi) - r(x))) < tol; 

% multiple roots
f = diskfun(@(t,r) cos(5*r), 'polar'); 
r = roots(f); 
rt = pi/5*[1, 2] - pi/10;
rt = rt.*exp(1i*x*pi); 

r1 = r(:,1); r2 = r(:,2); 
pass(4) = (norm(rt(:,1) - r1(x)))<tol; 
pass(5) = (norm(rt(:,2) - r2(x)))<tol; 

% pair of functions: 
x = diskfun(@(x,y) x); 
y = diskfun(@(x,y) y); 
f = x.^2 + y.^2 - 1/4; 
r = roots(f,x); 
pass(6) = (norm(r - [0, -.5; 0, .5]))<tol; 
end
