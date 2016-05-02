function pass = test_fevalm( pref ) 
% Test chebfun2/fevalm 

if ( nargin == 0) 
    pref = chebfunpref; 
end

rng(2016);
tol = 100*pref.cheb2Prefs.chebfun2eps;

% Check empty chebfun2: 
f = chebfun2; 
s = 2*rand(5,1) - 1; 
t = 2*rand(5,1) - 1; 
B = fevalm(f, s, t); 
pass(1) = isempty( B ); 

% Check symmetric function: 
f = chebfun2(@(x,y) cos(x.*y)); 
s = 2*rand(5,1) - 1; 
t = 2*rand(5,1) - 1; 
[ss, tt] = meshgrid( s, t); 
A = feval(f, ss, tt); 
B = fevalm(f, s, t); 
pass(2) = norm( A - B ) < tol; 

% Check essentially one dimensional function:
f = chebfun2(@(x,y) cos(y-.1)); 
s = 2*rand(5,1) - 1; 
t = 2*rand(5,1) - 1; 
[ss, tt] = meshgrid( s, t); 
A = feval(f, ss, tt); 
B = fevalm(f, s, t); 
pass(3) = norm( A - B ) < tol; 

% Check complex-valued function:
f = chebfun2(@(z) cos(z)); 
s = 2*rand(6,1) - 1; 
t = 2*rand(6,1) - 1; 
[ss, tt] = meshgrid( s, t); 
A = feval(f, ss, tt); 
B = fevalm(f, s, t); 
pass(4) = norm( A - B ) < tol; 

end