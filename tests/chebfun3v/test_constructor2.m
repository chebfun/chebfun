function pass = test_constructor2(pref)
% Test the Chebfun3v constructor when performing simple arithmetic
% operations.

if ( nargin < 1 )
    pref = chebfunpref;
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

% check the vectorize flag: 
f1 = @(x,y,z) x.*y.*z;
f2 = @(x,y,z) x*y*z;
H1 = chebfun3v(f1, f1);
H2 = chebfun3v(f2, f2, 'vectorize');
pass(1) = norm(H1 - H2) < tol; 

f1 = @(x,y,z) x.*y.*z; 
f2 = @(x,y,z) x*y*z; 
dom = [-2 3 -1 0 -1 1];
H1 = chebfun3v(f1, f1, dom);
H2 = chebfun3v(f2, f2, 'vectorize', dom);
pass(2) = norm(H1 - H2) < tol; 

f1 = @(x,y,z) x.*y.*z;
f2 = @(x,y,z) x*y*z;
H1 = chebfun3v(f1, f1, f1);
H2 = chebfun3v(f2, f2, f2, 'vectorize');
pass(3) = norm(H1 - H2) < tol; 

f1 = @(x,y,z) x.*y.*z;
f2 = @(x,y,z) x*y*z;
H1 = chebfun3v(f1, f1, f1, dom);
H2 = chebfun3v(f2, f2, f2, 'vectorize', dom);
pass(4) = norm(H1 - H2) < tol;

end