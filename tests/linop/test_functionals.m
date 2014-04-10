function pass = test_functionals
% Tests the construction and application of some functionals. 
% 
%  Toby Driscoll, 3 Feb 2014

d = [-1,2];
x = chebfun(@(x) x, d);
f = cos(x)./(1+x.^2);

%%
[z,e,s,dt] = linop.primitiveFunctionals(d);
A = [s-z; -2*dt(x.^2)];
Af = A*chebmatrix(f);
pass(1) = abs( sum(f) - Af{1} ) < 100*eps;
pass(2) = abs( -2*(f'*x.^2) - Af{2} ) < 100*eps;

%%
A = [e(2);e(0)];
Af = A*chebmatrix(f);
pass(3) = abs( f(2) - Af{1} ) < 100*eps;
pass(4) = abs( f(0) - Af{2} ) < 100*eps;

end
