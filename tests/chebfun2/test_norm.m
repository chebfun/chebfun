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

%% Check more norm() syntax. 
f = chebfun2(@(x,y) 1+x.^2.*y.^2 );
s = svd(f);
% Frobeius norm: 
pass(7) = abs( norm( f ) - sqrt(sum(s.^2)) ) < tol; 

% Frobeius norm:
pass(8) = abs( norm( f, 'fro' ) - sqrt(sum(s.^2)) ) < tol; 

% Operator norm, largest singular value: 
pass(9) = abs( norm( f, 'op' ) - s(1) ) < tol; 

% inf norm: 
pass(10) = abs( norm( f, 'inf' ) - 2 ) < tol; 
pass(11) =  abs( norm( f, inf ) - 2 ) < tol;

% max norm: 
pass(12) = abs( norm( f, 'max'  ) - 2 ) < tol;  
end