function pass = test_svd( ) 
% Test diskfun svd() command 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = diskfun(@(x,y) cos(x.*y.^2) ); 
s = svd( f ); 
pass(1) = abs( norm(f).^2 - sum(s.^2) ) < tol ; 

% Check resolved: 
pass(2) = ( s(end) < 1e1*tol );

% Scale invariant: 
g = 100*f; 
t = svd( g ); 
pass(3) = norm( s - t/100 ) < tol; 

end 