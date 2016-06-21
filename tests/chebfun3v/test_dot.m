function pass = test_dot(pref)
% Test DOT
if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 50*pref.cheb3Prefs.chebfun3eps;

% Check definition: 
F = chebfun3v(@(x,y,z) cos(x), @(x,y,z) sin(y), @(x,y,z) exp(z));
G = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
dotF1 = dot(F, G);
dotF2 = F' * G;
pass(1) = norm(dotF1 - dotF2) < tol;

% Check definition again:
F = chebfun3v(@(x,y,z) cos(x), @(x,y,z) sin(y), @(x,y,z) x.*y);
G = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
dotF1 = dot(F, G);
dotF2 = F' * G;
pass(2) = norm(dotF1 - dotF2) < tol;

end