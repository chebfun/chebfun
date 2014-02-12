function pass = test_Coefficients( pref )
% Test to check that Chebfun2 can compute its bivariate tensor Chebyshev 
% coefficients correctly.

if ( nargin < 1 ) 
    pref = chebpref; 
end 
tol = 100 * pref.cheb2Prefs.eps; 

n = 10; 
f = chebpoly(n); 
g = chebpoly(n); 
h = chebfun2(@(x,y) f(x).*g(y)); 
X = chebpoly2(h); 

pass(1)  = (abs( X(end-n,end-n) -1) < tol );  
X(end-n,end-n) = X(end-n,end-n)-1;
pass(2) =  (norm(X) < tol ); 

% Are cheb2poly and cheb2polyval inverses?  
f=chebfun2( @(x,y) cos(x+y) ); 
lenc = length(f.cols); 
lenr = length(f.rows);
[xx, yy] = meshgrid(chebpts(lenc),chebpts(lenr)); 
vals = f(xx, yy); 
X = chebfun2.coeffs2vals( chebpoly2( f ) ); 
pass(3) = ( (norm(X-vals)) < tol );

end