function pass = test_adaptivity( prefs )
% Test to check adaptivity decisions can be made by the user. 

if ( nargin < 1 ) 
    prefs = chebfunpref(); 
end 
tol = 100*prefs.cheb2Prefs.eps; 


d = [-5 1 0 .1];
N = chebop2(@(u) diffy(u) + diffx(u,3),d);
N.dbc = @(x) exp(-10*x.^2);
N.rbc = @(t,u) [u diff(u)];
N.lbc = 0;
u1 = mldivide(N, 0, 200, inf);
[m, n] = length( u1 ); 
pass(1) = ( n == 200 ); 

N = chebop2(@(u) diffy(u) + diffx(u,3),d);
N.dbc = @(x) exp(-10*x.^2);
N.rbc = @(t,u) [u diff(u)];
N.lbc = 0;
u2 = mldivide(N, 0, 200, 200);
[m, n] = length( u2 ); 
pass(2) = ( m == 200 ); 
pass(3) = ( n == 200 ); 

% check solutions are about the same (PDE is not very numerically stable 
% because bcs do not match): 
pass(4) = norm( u1 - u2 ) < 20*sqrt(tol); 

end