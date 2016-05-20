function pass = test_conj(pref)
% Test CONJ

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

% Test whether the conjugate of a real-valued CHEBFUN3 is itself.
f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
h = conj(f);
pass(1) = norm(f - h) < tol;

% Test whether the conjugate of 1i*f is -1i*f.
f = chebfun3(@(x,y,z) cos(x.*y.*z));
h = conj(1i*f);
pass(2) = norm(1i*f + h) < 100*tol;

% Test whether the conjugate of f+1i*g is f-1i*g.
f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
g = chebfun3(@(x,y,z) sin(x + y.^2 + z.^3));
h = conj(f + 1i*g); 
pass(3) = norm(f - 1i*g - h) < tol;

end