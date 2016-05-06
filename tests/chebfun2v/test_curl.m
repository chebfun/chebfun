function pass = test_curl( pref ) 
% Test CURL
if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 50*pref.cheb2Prefs.chebfun2eps;

% Check definition: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y));  
curlF = diff(F(2), 1, 2) - diff(F(1), 1, 1);
pass(1) = ( norm(curlF - curl(F)) < tol );

% Check definition: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y), @(x,y) x.*y); 
curlF = [ diff(F(3),1,1) ; -diff(F(3),1,2) ;...
          diff(F(2),1,2) - diff(F(1),1,1) ];
pass(2) = ( norm(curlF - curl(F)) < tol );

% curl of a gradient field is zero: 
f = chebfun2(@(x,y) cos(x+y.^2) + sin(y) + y); 
F = gradient( f ); 
pass(3) = ( norm( curl( F ) ) < tol ); 

end
