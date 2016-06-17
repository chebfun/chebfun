function pass = test_subsref(pref)
% Test Chebfun3v subsref() command.

if ( nargin == 0) 
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

% Test recursive subsref:  
f = chebfun3(@(x,y,z) sin(x.*y.*z));
F = [f; f; f];

G = F(1);
exactCore = G.core; 

pass(1) = norm(exactCore(:) - F(1).core(:)) < tol;

end