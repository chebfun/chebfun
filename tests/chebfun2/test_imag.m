function pass = test_imag( pref ) 
% Test IMAG
if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 100*pref.eps; 

f = chebfun2(@(x,y) cos(x.*y)); 
g = chebfun2(@(x,y) sin(x+y.^2));

% Simple consistency check: 
x = linspace(-1,1); 
[xx, yy] = meshgrid(x);
h = f + 1i *g; 
pass(1) = norm( imag( feval(h,xx,yy) ) - g(xx,yy) ) < 10*tol;
pass(2) = norm( imag( h ) - g ) < tol;

end