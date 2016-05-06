function pass = test_norm( ) 
% Test for the norm of a chebfun2: 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = chebfun2(@(x,y) 1+x.^2.*y.^2 );
s = svd(f);
% Frobeius norm: 
pass(1) = abs( norm( f ) - sqrt(sum(s.^2)) ) < tol; 

% Frobeius norm:
pass(2) = abs( norm( f, 'fro' ) - sqrt(sum(s.^2)) ) < tol; 

% Operator norm, largest singular value: 
pass(3) = abs( norm( f, 'op' ) - s(1) ) < tol; 

% inf norm: 
pass(4) = abs( norm( f, 'inf' ) - 2 ) < tol; 
pass(5) =  abs( norm( f, inf ) - 2 ) < tol;

% max norm: 
pass(5) = abs( norm( f, 'max'  ) - 2 ) < tol;  

end