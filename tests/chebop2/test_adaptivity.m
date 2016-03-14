function pass = test_adaptivity( pref )
% Test to check adaptivity decisions can be made by the user. 

if ( nargin < 1 ) 
    pref = chebfunpref();
end
tol = 100*pref.cheb2Prefs.chebfun2eps;

% This is the first test so remove the warning: 
state = warning; 
warning('off','CHEBFUN:CHEBOP2:chebop2:experimental');

d = [-5 1 0 .1];
N = chebop2(@(u) diffy(u) + diffx(u,3),d);
N.dbc = @(x) exp(-10*x.^2);
N.rbc = @(t,u) [u ; diff(u)];
N.lbc = 0;
u1 = mldivide(N, 0, 200, inf);
[m, n] = length( u1 ); 
pass(1) = ( m == 200 ); 

u2 = mldivide(N, 0, 200, 200);
[m, n] = length( u2 ); 
pass(2) = ( m == 200 ); 
pass(3) = ( n == 200 ); 

u3 = mldivide(N, 0, inf, 200);
[m, n] = length( u2 ); 
pass(4) = ( n == 200 ); 

% check small values of n: 
u4 = mldivide(N, 0, 10, 10);
[m, n] = length( u4 ); 
pass(5) = ( m == 10 ); 
pass(6) = ( n == 10 );  

% check small values of n: 
u4 = mldivide(N, 0, 6, 10);
[m, n] = length( u4 ); 
pass(7) = ( n == 10 ); 
pass(8) = ( m == 6 );

% check solutions are about the same (PDE is not very numerically stable 
% because bcs do not match): 
pass(9) = norm( u1 - u2 ) < 20*sqrt(tol); 
pass(10) = norm( u2 - u3 ) < 20*sqrt(tol); 

% restore the warnings: 
warning(state);

end
