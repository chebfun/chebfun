function pass = test_integral3(pref)
% Test definite integral of f over a restricted cuboid.

if ( nargin == 0)
    pref = chebfunpref; 
end
tol = 1e4*pref.cheb3Prefs.chebfun3eps;

% Check empty case
pass(1) = isempty(integral3(chebfun3()));

% Example: 
f = chebfun3(@(x,y,z) x, [0, 2, 0, 3, -1, 4]);
I = integral3(f, [0, 1, 0, 2, 0, 3]);
exact = 3;
pass(2) = abs(I - exact) < tol;

end