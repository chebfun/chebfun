function pass = test_fevalm( pref ) 
% Test diskfun/fevalm 

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 100*pref.cheb2Prefs.chebfun2eps;
rng(2016);

% Check empty diskfun: 
f = diskfun;
t = pi*(2*rand(5,1) - 1); 
r = 2*rand(5,1)-1; 
B = fevalm(f, t, r); 
pass(1) = isempty( B );

% Check rank 1 diskfun: 
f = diskfun(@(t,r) r.*sin(t), 'polar'); 
t = pi*(2*rand(5,1) - 1); 
r = rand(5,1); 
[tt, rr] = meshgrid( t, r); 
A = feval(f, tt, rr, 'polar'); 
B = fevalm(f,t,r); 
pass(2) = norm( A - B ) < tol; 

% Check essentially one dimensional function:
f = diskfun(@(t,r) exp(-r.^2), 'polar'); 
t = pi*(2*rand(5,1) - 1); 
r = rand(5,1); 
[tt, rr] = meshgrid( t, r); 
A = feval(f, tt, rr, 'polar'); 
B = fevalm(f, t, r); 
pass(3) = norm( A - B ) < tol; 


end