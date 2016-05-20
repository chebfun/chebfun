function pass = test_curl(pref)
% Test CURL
if ( nargin == 0 ) 
    pref = chebfunpref; 
end
tol = 50*pref.cheb3Prefs.chebfun3eps;

% Check definition: 
F = chebfun3v(@(x,y,z) cos(x), @(x,y,z) sin(y), @(x,y,z) x.*y); 
curlF = [ diff(F(3), 1, 2) - diff(F(2), 1, 3);
          diff(F(1), 1, 3) - diff(F(3), 1, 1);
          diff(F(2), 1, 1) - diff(F(1), 1, 2) ];
pass(1) = norm(curlF - curl(F)) < tol;

% Curl of a gradient vector field is zero:
f = chebfun3(@(x,y,z) cos(x+y.^2) + sin(y) + y); 
F = gradient(f);
pass(2) = norm(curl(F)) < tol;

end