function pass = test_zerofunction( pref )
% This test checks that a zero chebfun2 is being treated correctly for
% common commands. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 100 * pref.cheb2Prefs.chebfun2eps;

% construction
f = chebfun2(0); 
f = chebfun2(@(x) 0*x); 
f = chebfun2(0,[-2 2 -3 3]); 

% adding together a 0 chebfun2. 
f = chebfun2(0); 
g = chebfun2(@(x,y) cos(x.*y)); 
pass(1) = (abs(norm(g + f) - norm(g)) < tol); 

v = abs(f(pi/6,pi/6)); 
v = v+abs(length(f)); 
v = v+abs(rank(f));
pass(2) = (v==1); 

v = abs(norm(sum(f)));
v = v + abs(integral2(f));

v = v + abs(norm(diff(f)));

v = v + abs(norm(diff(f,1,1)));

v = v + abs(svd(f));

pass(3) = (v==0); 


% evaluation on an array. 
r = rand(10,8); 
v = f(r,r); 
if all(size(v) == size(r))
    pass(4) = 1; 
end

f = chebfun2( zeros(5,4) ); 
pass(5) = ( norm( coeffs2(f) - zeros(5,4) ) == 0 );

f = chebfun2( 0 ); 
[C, D, R] = coeffs2( f, 11, 13 );
pass(6) = ( norm( C - zeros(11,1) ) == 0 && ...
            norm( D - 1 )           == 0 && ...
            norm( R - zeros(13,1) ) == 0 );

end