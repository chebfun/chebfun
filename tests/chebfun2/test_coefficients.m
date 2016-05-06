function pass = test_coefficients( pref )
% Test to check that Chebfun2 can compute its bivariate tensor Chebyshev 
% coefficients correctly.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 100 * pref.cheb2Prefs.chebfun2eps;

kind = 2;
if ( isa( pref.tech(), 'chebtech2' ) )
    kind = 2; 
elseif ( isa( pref.tech(), 'chebtech1' ) )
    kind = 1; 
end

n = 10; 
f = chebpoly(n); 
g = chebpoly(n); 
h = chebfun2(@(x,y) f(x).*g(y)); 
X = chebcoeffs2(h); 

pass(1)  = (abs( X(n+1,n+1) - 1) < tol );  
X(n+1,n+1) = X(n+1,n+1)-1;
pass(2) =  (norm(X) < tol ); 

% Are cheb2poly and cheb2polyval inverses?  
f = chebfun2( @(x,y) cos(x+y) ); 
lenc = length(f.cols); 
lenr = length(f.rows);
[xx, yy] = chebfun2.chebpts2(lenc,lenr,[-1,1,-1,1],kind); 
vals = f(xx, yy); 
X = chebfun2.coeffs2vals( chebcoeffs2( f ) ); 
pass(3) = ( (norm(X-vals)) < tol );

end