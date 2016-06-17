function pass = test_integral2(pref)
% Test surface integral for 3D functions.

if ( nargin == 0)
    pref = chebfunpref; 
end
tol = 1e7*pref.cheb3Prefs.chebfun3eps;

% Example: 
S = chebfun2v(@(u,v) u.*cos(v), @(u,v) u.*sin(v), @(u,v) v, [0 4 0 2*pi]);
f = chebfun3(@(x,y,z) sqrt(1+x.^2 + y.^2), [-4 4 -4 4 0 2*pi]);
I = integral2(f, S);
exact = 152*pi/3;       % From https://www.youtube.com/watch?v=pAWLCFYsrVs
pass(1) = abs(I - exact) < tol;

end