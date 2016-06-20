function pass = test_integral2(pref)
% Test chebfun3v/integral2

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1e4*pref.cheb3Prefs.chebfun3eps;

% Surface parametrisation (unit disk in xy-plane)
S = chebfun2v(@(r,phi) r.*cos(phi), @(r,phi) r.*sin(phi), @(r,phi) 0*r, ... 
    [0, 1, 0, 2*pi]);

% Constant vector field [0; 0; 1]
F = chebfun3v(@(x,y,z) 0*x, @(x,y,z) 0*y, @(x,y,z) 0*x + 1);

% Flux integral, exact value is pi:
I = integral2(F, S);

pass(1) = ( abs(I-pi) < tol );

end