function pass = test_integral(pref)
% Test path integral of 3D vector fields:

if ( nargin == 0)
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb3Prefs.chebfun3eps;

[x,y,z] = cheb.xyz;

% http://tutorial.math.lamar.edu/Classes/CalcIII/LineIntegralsVectorFields.aspx
F = [8*x.^2.*y.*z; 5*z; -4*x.*y];
curve = chebfun(@(t) [t, t.^2, t.^3], [0, 1]);
exact = 1;
pass(1) = abs(integral(F, curve) - exact) < tol;

% http://tutorial.math.lamar.edu/Classes/CalcIII/LineIntegralsPtII.aspx
dom = [-1, 1, -1, 1, 0, 4*pi^2];
x = chebfun3(@(x,y,z) x, dom);
y = chebfun3(@(x,y,z) y, dom);
z = chebfun3(@(x,y,z) z, dom);
F = [y; x; z];
curve = chebfun(@(t) [cos(t), sin(t), t.^2], [0, 2*pi]);
exact = 8*pi^4;
pass(2) = abs(integral(F, curve) - exact)/exact < tol;

end