function pass = test_coefficients(pref)
% Test whether Chebfun3 can compute its bivariate tensor Chebyshev 
% coefficients correctly.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 100 * pref.cheb3Prefs.chebfun3eps;

kind = 2;
if ( isa(pref.tech(), 'chebtech2') )
    kind = 2;
elseif ( isa( pref.tech(), 'chebtech1') )
    kind = 1; 
end

n = 10; 
f = chebpoly(n);
g = chebpoly(n);
h = chebpoly(n);
k = chebfun3(@(x,y,z) f(x).*g(y).*h(z));
X = chebcoeffs3(k);

pass(1) = abs( X(n+1, n+1, n+1) - 1) < tol;
X(n+1, n+1, n+1) = X(n+1, n+1, n+1) - 1;
pass(2) =  norm(X(:)) < tol;

% Are cheb2poly and cheb2polyval inverses?  
f = chebfun3( @(x,y,z) cos(x+y+z) );
lenCols = length(f.cols);
lenRows = length(f.rows);
lenTubes = length(f.tubes);
[xx, yy, zz] = chebpts3(lenCols, lenRows, lenTubes, ...
    [-1, 1, -1, 1, -1, 1], kind);
vals = f(xx, yy, zz);
X = chebfun3.coeffs2vals(chebcoeffs3(f));
pass(3) = (norm(X(:)-vals(:))) < tol;

end