function pass = test_chebfun2v_ctor( pref ) 
% Test the Chebfun2v constructor when performing simple arithmetic
% operations. 

if ( nargin < 1 ) 
    pref = chebpref; 
end 

pass = 1; 

try 
% % Adaptive calls % % 
f = @(x,y) cos(x); 
f = chebfun2v(f, f); 

g = @(x,y) sin(y); 
g = chebfun2v(g, g); 

% exact answers. 
plus_exact = @(x,y) cos(x) + sin(y); 
chebfun2v(plus_exact, plus_exact); 

minus_exact = @(x,y) cos(x) - sin(y); 
chebfun2v(minus_exact, minus_exact); 

mult_exact = @(x,y) cos(x).*sin(y); 
chebfun2v(mult_exact, mult_exact); 

pow_exact = @(x,y) cos(x).^sin(y); 
chebfun2v(pow_exact, pow_exact); 

catch
    pass = 0 ; 
end
end