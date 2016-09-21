function pass = test_gradient(pref)
% Check the gradient command in Chebfun3.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e4 * pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) x);
exactGrad = chebfun3v(@(x,y,z) 1+0*x, @(x,y,z) 0*x, @(x,y,z) 0*x);
fGrad = grad(f);
pass(1) = norm(exactGrad - fGrad) < tol; 

f = chebfun3(@(x,y,z) cos(x).*exp(y).*sin(z));
exactGrad = chebfun3v(@(x,y,z) -sin(x).*exp(y).*sin(z), ...
    @(x,y,z) cos(x).*exp(y).*sin(z), @(x,y,z) cos(x).*exp(y).*cos(z));
fGrad = grad(f);
pass(2) = norm(exactGrad - fGrad) < tol; 

f = chebfun3(@(x,y,z) cos(x.*y.*z));
exactGrad = chebfun3v(@(x,y,z) -y.*z.*sin(x.*y.*z), ...
    @(x,y,z) -x.*z.*sin(x.*y.*z), @(x,y,z) -x.*y.*sin(x.*y.*z));
fGrad = grad(f);
pass(3) = norm(exactGrad - fGrad) < tol; 

f = chebfun3(@(x,y,z) x.^2 + x.*y.^2.*z.^3);
exactGrad = chebfun3v(@(x,y,z) 2*x+y.^2.*z.^3, @(x,y,z) 2*x.*y.*z.^3, ...
    @(x,y,z) 3*x.*y.^2.*z.^2);
fGrad = grad(f);
pass(4) = norm(exactGrad - fGrad) < tol; 

end