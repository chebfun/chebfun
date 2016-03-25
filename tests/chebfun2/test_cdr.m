function pass = test_cdr( pref ) 
% Test CDR 

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb2Prefs.chebfun2eps;
j = 1; 

f = chebfun2(@(x,y) cos(x.*y));
x = linspace(-1,1)';
[xx,yy] = meshgrid(x);
[C, D, R] = cdr( f ); 
err = norm( f(xx, yy) - C(x,:)*D*R(x,:).' );
pass(j) = err < tol; j = j + 1; 

f = chebfun2(@(x,y) cos(x.*y), [-3 4 -1 10]);
x = linspace(-3,4)';
y = linspace(-1,10)';
[xx,yy] = meshgrid(x, y);
[C, D, R] = cdr( f ); 
err = norm( f(xx, yy) - C(y,:)*D*R(x,:).' );
pass(j) = err < tol; j = j + 1; 

d = cdr( f ); 
pass(j) = norm( diag( D ) - d(:) ) < tol; j = j + 1; 

end