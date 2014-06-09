function pass = test_lagpoly( prefs ) 
% Test the LAGPOLY command. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end
tol = 1e5*prefs.techPrefs.eps;

% Laguerre polynomials:
h0 = lagpoly(0); 
pass(1) = norm( h0 - chebfun(@(x) 1+0*x, [0,inf] ) ) < tol;

h1 = lagpoly(1); 
pass(2) = norm( h1 - chebfun(@(x) -x+1, [0,inf] ) ) < tol;

h2 = lagpoly(2); 
x = linspace(0,1);
pass(3) = norm( feval(h2,x) - (x.^2-4*x+2)/2 ) < tol;

h3 = lagpoly(3); 
pass(4) = norm( feval(h3,x) - (-x.^3+9*x.^2-18*x+6)/6 ) < tol;

h4 = lagpoly(4); 
pass(5) = norm( feval(h4,x) - (x.^4-16*x.^3+72*x.^2-96*x+24)/24 ) < tol;

h5 = lagpoly(5); 
pass(6) = norm( feval(h5,x) -...
                  (-x.^5+25*x.^4-200*x.^3+600*x.^2-600*x+120)/120 ) < 10*tol;

h6 = lagpoly(6); 
pass(7) = norm( feval(h6,x) -...
       (x.^6-36*x.^5+450*x.^4-2400*x.^3+5400*x.^2-4320*x+720)/720 ) < 100*tol;