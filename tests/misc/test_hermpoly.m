function pass = test_hermpoly( prefs ) 
% Test the HERMPOLY command. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end
tol = 1e5*prefs.techPrefs.eps;

% Physicists' Hermite polynomials
h0 = hermpoly(0); 
pass(1) = norm( h0 - chebfun(@(x) 1+0*x, [-inf,inf] ) ) < tol;

h1 = hermpoly(1); 
pass(2) = norm( h1 - chebfun(@(x) 2*x, [-inf,inf] ) ) < tol;

h2 = hermpoly(2); 
pass(3) = norm( h2 - chebfun(@(x) 4*x.^2-2, [-inf,inf] ) ) < tol;

h3 = hermpoly(3); 
x = linspace(-1,1);
pass(4) = norm( feval(h3,x) - 8*x.^3+12*x ) < tol;

h4 = hermpoly(4); 
x = linspace(-1,1);
pass(5) = norm( feval(h4,x) - (16*x.^4-48*x.^2+12) ) < tol;

h5 = hermpoly(5); 
x = linspace(-1,1);
pass(6) = norm( feval(h5,x) - (32*x.^5-160*x.^3+120*x) ) < 10*tol;

h6 = hermpoly(6); 
x = linspace(-1,1);
pass(7) = norm( feval(h6,x) - (64*x.^6-480*x.^4+720*x.^2-120) ) < 100*tol;


% probabilists' Hermite polynomials:
h0 = hermpoly(0, 'prob'); 
pass(8) = norm( h0 - chebfun(@(x) 1+0*x, [-inf,inf] ) ) < tol;

h1 = hermpoly(1, 'prob'); 
pass(9) = norm( h1 - chebfun(@(x) x, [-inf,inf] ) ) < tol;

h2 = hermpoly(2, 'prob'); 
pass(10) = norm( h2 - chebfun(@(x) x.^2-1, [-inf,inf] ) ) < tol;

h3 = hermpoly(3, 'prob'); 
x = linspace(-1,1);
pass(11) = norm( feval(h3,x) - (x.^3-3*x) ) < tol;

h4 = hermpoly(4, 'prob'); 
x = linspace(-1,1);
pass(12) = norm( feval(h4,x) - (x.^4-6*x.^2+3) ) < tol;

h5 = hermpoly(5, 'prob'); 
x = linspace(-1,1);
pass(13) = norm( feval(h5,x) - (x.^5-10*x.^3+15*x) ) < tol;

h6 = hermpoly(6, 'prob'); 
x = linspace(-1,1);
pass(14) = norm( feval(h6,x) - (x.^6-15*x.^4+45*x.^2-15) ) < 10*tol;
