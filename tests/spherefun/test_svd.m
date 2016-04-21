function pass = test_svd( ) 
% Test spherefun svd() command 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = spherefun(@(x,y,z) cos(x.*y.*z) ); 
s = svd( f ); 
pass(1) = abs( norm(f).^2 - sum(s.^2) ) < tol ; 

% Check resolved: 
pass(2) = ( s(end) < tol );

% Scale invariant: 
g = 100*f; 
t = svd( g ); 
pass(3) = norm( s - t/100 ) < tol; 

end 