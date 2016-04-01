function pass = test_BMCsvd( ) 
% Test spherefun BMCsvd() command 

tol = 1e3*chebfunpref().cheb2Prefs.chebfun2eps;

f = spherefun(@(x,y,z) sin(pi*x.*y) + sin(pi*x.*z));
s = BMCsvd( f ); 

% Check resolved: 
pass(1) = ( s(end) < 1e3*tol );

% Scale invariant: 
g = 100*f; 
t = BMCsvd( g ); 
pass(2) = norm( s - t/100 ) < tol; 

end 