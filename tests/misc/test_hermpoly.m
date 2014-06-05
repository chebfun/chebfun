function pass = test_hermpoly( prefs ) 
% Test the hermpoly command. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end
tol = 1e5*prefs.techPrefs.eps;

h0 = hermpoly(0); 
pass(1) = norm( h0 - chebfun(@(x) 1+0*x, [-inf,inf] ) ) < tol;

h1 = hermpoly(1); 
pass(1) = norm( h1 - chebfun(@(x) 2*x, [-inf,inf] ) ) < tol;

h2 = hermpoly(2); 
pass(1) = norm( h2 - chebfun(@(x) 4*x.^2-2, [-inf,inf] ) ) < tol;

h3 = hermpoly(3); 
x = linspace(-1,1);
pass(1) = norm( feval(h3,x) - 8*x.^3+12*x ) < tol;

h4 = hermpoly(4); 
x = linspace(-1,1);
pass(1) = norm( feval(h4,x) - (16*x.^4-48*x.^2+12) ) < tol;

h5 = hermpoly(5); 
x = linspace(-1,1);
pass(1) = norm( feval(h5,x) - (32*x.^5-160*x.^3+120*x) ) < tol;

h6 = hermpoly(6); 
x = linspace(-1,1);
pass(1) = norm( feval(h6,x) - (64*x.^6-480*x.^4+720*x.^2-120) ) < tol;
