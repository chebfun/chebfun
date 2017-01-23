function pass = test_norm(pref)

if ( nargin < 1 )
    pref = chebfunpref; 
end 
tol = 1000*pref.cheb2Prefs.chebfun2eps;

%% Real-valued function:
f = chebfun2(@(x,y) x);
exact = sqrt(4/3);      % sqrt(sum2(x.^2))
pass(1) = abs(norm(f) - exact) < tol;

exact = 1;
pass(2) = abs(norm(f, 'inf') - exact) < tol;

% p-norm for even values of p should work
exact = (4/5)^(1/4);
pass(3) = abs(norm(f, 4) - exact) < tol;

%% Complex-valued function:
f = chebfun2(@(x,y) 1i*x);
exact = sqrt(4/3);      % sqrt(sum2(x.^2))
pass(4) = abs(norm(f) - exact) < tol;

exact = 1;
pass(5) = abs(norm(f, 'inf') - exact) < tol;

% p-norm for even values of p should work
exact = (4/5)^(1/4);
pass(6) = abs(norm(f, 4) - exact) < tol;

%% functions with rank>1
f = chebfun2(@(x,y) exp(x.*y));

exact = 2.236768845167052;
pass(7) = abs(norm(f) - exact) < tol; % default is the Frobenius norm 
pass(8) = abs(norm(f,'fro') - exact) < tol; % Frobenius norm 

exact =  2.119814813637055;
pass(7) = abs(norm(f,2) - exact) < tol; % spectral (operator) norm 

exact =  2.925303491814361;
pass(8) = abs(norm(f,'nuc') - exact) < tol; % Frobenius norm 

end