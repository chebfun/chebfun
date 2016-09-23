function pass = test_divergence(pref)
% Test DIVERGENCE

if ( nargin == 0 ) 
    pref = chebfunpref;
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

% Check definition:
F = chebfun3v(@(x,y,z) cos(x), @(x,y,z) sin(y), @(x,y,z) exp(z));
divF = diff(F(1), 1, 1) + diff(F(2), 1, 2) + diff(F(3), 1, 3);
pass(1) = norm(divF - divergence(F)) < tol;

% Check divergence of gradient is laplacian: 
f = chebfun3(@(x,y,z) cos(x.*y.*z));
pass(2) = norm(laplacian(f) - divergence(gradient(f))) < tol;

end