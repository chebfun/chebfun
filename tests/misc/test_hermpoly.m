function pass = test_hermpoly( prefs ) 
% Test the HERMPOLY command. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end
tol = 1e6*prefs.techPrefs.chebfuneps;

% Physicists' Hermite polynomials
h0 = hermpoly(0); 
err = norm( h0 - chebfun(@(x) 1+0*x, [-inf,inf] ) );
pass(1) = err < tol;

h1 = hermpoly(1); 
err = norm( h1 - chebfun(@(x) 2*x, [-inf,inf] ) );
pass(2) = err < tol;

h2 = hermpoly(2); 
err = norm( h2 - chebfun(@(x) 4*x.^2-2, [-inf,inf] ) );
pass(3) = err < tol;

h3 = hermpoly(3); 
x = linspace(-1,1);
err = norm( feval(h3,x) - 8*x.^3+12*x );
pass(4) = err < tol;

h4 = hermpoly(4); 
x = linspace(-1,1);
err = norm( feval(h4,x) - (16*x.^4-48*x.^2+12) );
pass(5) = err < tol;

h5 = hermpoly(5); 
x = linspace(-1,1);
err = norm( feval(h5,x) - (32*x.^5-160*x.^3+120*x) );
pass(6) = err < 10*tol;

h6 = hermpoly(6); 
x = linspace(-1,1);
err = norm( feval(h6,x) - (64*x.^6-480*x.^4+720*x.^2-120) );
pass(7) = err < 100*tol;

% probabilists' Hermite polynomials:
h0 = hermpoly(0, 'prob'); 
err = norm( h0 - chebfun(@(x) 1+0*x, [-inf,inf] ) );
pass(8) = err < tol;

h1 = hermpoly(1, 'prob'); 
err = norm( h1 - chebfun(@(x) x, [-inf,inf] ) );
pass(9) = err < tol;

h2 = hermpoly(2, 'prob');
g2 = chebfun(@(x) x.^2-1, [-inf,inf] );
err = norm( h2 - chebfun(@(x) x.^2-1, [-inf,inf] ) );
pass(10) = err < tol;

h3 = hermpoly(3, 'prob'); 
x = linspace(-1,1);
err = norm( feval(h3,x) - (x.^3-3*x) );
pass(11) = err < tol;

h4 = hermpoly(4, 'prob'); 
x = linspace(-1,1);
err = norm( feval(h4,x) - (x.^4-6*x.^2+3) );
pass(12) = err < tol;

h5 = hermpoly(5, 'prob'); 
x = linspace(-1,1);
err = norm( feval(h5,x) - (x.^5-10*x.^3+15*x) );
pass(13) = err < tol;

h6 = hermpoly(6, 'prob'); 
x = linspace(-1,1);
norm( feval(h6,x) - (x.^6-15*x.^4+45*x.^2-15) );
pass(14) = err < 10*tol;
