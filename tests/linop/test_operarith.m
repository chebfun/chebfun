function pass = linop_operarith
% This test checks basic arithmetric operations of linops
% Toby Driscoll, 3 Feb 2014

d = [-1,4];
[Z,I,D,C,M] = linop.primitiveOperators(d);
f = chebfun(@(x) exp(sin(x).^2+2),d);

%%
F = M(f);
A = -(2*D^2 - F*C + 3*I);
Af = A*f;

%%
pass = norm( Af - (f.*cumsum(f)-2*diff(f,2)-3*f) ) < 1e4*eps;

end
