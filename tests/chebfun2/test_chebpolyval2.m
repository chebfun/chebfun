function pass = test_chebpolyval2()
% Check the chebpolyval2 commands in trunk and @chebfun2 folder 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 100 * pref.cheb2Prefs.chebfun2eps;

kind = 2;
if ( isa( pref.tech(), 'chebtech2' ) )
    kind = 2; 
elseif ( isa( pref.tech(), 'chebtech1' ) )
    kind = 1; 
end

% check the trunk chebpolyval2 command.
n = 20;
T = chebpoly(n);
N = 5*n;
[xx, yy] = chebfun2.chebpts2(N,N,[-1 1 -1 1],kind); 
A = T(xx).*T(yy);  

C = zeros(N); C(n+1,n+1)=1;   
X = chebfun2.coeffs2vals(C); 
pass(1) = ( norm(A - X) < 100*tol); 

% check the @chebfun2/chebpolyval2 command.
f = chebfun2(@(x,y) cos(x.*y) + sin(x) + exp(y) ); 
[A1, A2, A3] = chebpolyval2( f ); 
X = chebpolyval2( f ); 
pass(2) = ( norm(X - A1*A2*A3.') < tol);

% Check when degree in x \neq degree in y: 
A = rand(3,4);
f = chebfun2(A);
B = chebpolyval2(f);
pass(3) = ( norm(A - B) < tol);
end