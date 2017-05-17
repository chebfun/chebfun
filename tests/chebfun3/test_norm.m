function pass = test_norm(pref)

if ( nargin < 1 )
    pref = chebfunpref; 
end 
tol = 1000*pref.cheb3Prefs.chebfun3eps;

%% Real-valued function:
f = chebfun3(@(x,y,z) x);
exact = 2*sqrt(2/3);      % sqrt(sum3(x.^2))
pass(1) = abs(norm(f) - exact) < tol;

exact = 1;
pass(2) = abs(norm(f, 'inf') - exact) < tol;

% p-norm for even values of p should work
exact = (8/5)^(1/4);
pass(3) = abs(norm(f, 4) - exact) < tol;

%% Complex-valued function:
f = chebfun3(@(x,y,z) 1i*x);
exact = 2*sqrt(2/3);      % sqrt(sum3(x.^2))
pass(4) = abs(norm(f) - exact) < tol;

exact = 1;
pass(5) = abs(norm(f, 'inf') - exact) < tol;

% p-norm for even values of p should work
exact = (8/5)^(1/4);
pass(6) = abs(norm(f, 4) - exact) < tol;


end