function pass = test_fevalm( pref ) 
% Test spherefun/fevalm 

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 100*pref.cheb2Prefs.chebfun2eps;
rng(2016);

% Check empty spherefun: 
f = spherefun;
s = pi*(2*rand(5,1) - 1); 
t = pi/2*rand(5,1); 
B = fevalm(f, s, t); 
pass(1) = isempty( B );

% Check rank 1 spherefun: 
f = chebfun2(@(lam,th) cos(lam).*sin(th)); 
s = pi*(2*rand(5,1) - 1); 
t = pi/2*rand(5,1); 
[ss, tt] = meshgrid( s, t); 
A = feval(f, ss, tt); 
B = fevalm(f, s, t); 
pass(2) = norm( A - B ) < tol; 

% Check essentially one dimensional function:
f = spherefun(@(lam,th) exp(-(cos(th)-1).^2)); 
s = pi*(2*rand(5,1) - 1); 
t = pi/2*rand(5,1); 
[ss, tt] = meshgrid( s, t); 
A = feval(f, ss, tt); 
B = fevalm(f, s, t); 
pass(3) = norm( A - B ) < tol; 

end