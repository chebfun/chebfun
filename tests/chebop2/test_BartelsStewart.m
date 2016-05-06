function pass = test_BartelsStewart( pref )
% Test generalized Sylvester matrix equation solver

if ( nargin < 1 )
    pref = chebfunpref();
end
tol = 1e4*pref.cheb2Prefs.chebfun2eps;

n = 10; 
rng(0)
A = rand(n); 
B = rand(n); 
C = rand(n); 
D = rand(n); 
X = rand(n); 

E = A * X * B.' + C * X * D.';

Y = chebop2.bartelsStewart(A, B, C, D, E, 0, 0); 
pass(1) = norm( Y - X ) < tol; 


A = rand(n) + 1i*rand(n); 
B = rand(n) + 1i*rand(n); 
C = rand(n) + 1i*rand(n); 
D = rand(n) + 1i*rand(n); 
X = rand(n) + 1i*rand(n); 

E = A * X * B.' + C * X * D.';

Y = chebop2.bartelsStewart(A, B, C, D, E, 0, 0); 
pass(2) = norm( Y - X ) < 10*tol; 


tol = 1000*tol; 
n = 100; 
rng(0)
A = rand(n); 
B = rand(n); 
C = rand(n); 
D = rand(n); 
X = rand(n); 

E = A * X * B.' + C * X * D.';

Y = chebop2.bartelsStewart(A, B, C, D, E, 0, 0); 
pass(3) = norm( Y - X ) < 10*tol; 


A = rand(n) + 1i*rand(n); 
B = rand(n) + 1i*rand(n); 
C = rand(n) + 1i*rand(n); 
D = rand(n) + 1i*rand(n); 
X = rand(n) + 1i*rand(n); 

E = A * X * B.' + C * X * D.';

Y = chebop2.bartelsStewart(A, B, C, D, E, 0, 0); 
pass(4) = norm( Y - X ) < 10*tol; 

end