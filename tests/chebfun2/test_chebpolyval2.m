function pass = test_chebpolyval2()
% Check the chebpolyval2 commands in trunk and @chebfun2 folder 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 100 * pref.eps; 

kind = 2;
if ( isa( pref.tech(), 'chebtech2' ) )
    kind = 2; 
elseif ( isa( pref.tech(), 'chebtech1' ) )
    kind = 1; 
end

% check the trunk chebpolyval2 command.
T = chebpoly(20); 
[xx, yy] = chebfun2.chebpts2(100,100,[-1 1 -1 1],kind); 
A = T(xx).*T(yy);  

C = zeros(100); C(end-20,end-20)=1;   
X = chebfun2.coeffs2vals(C); 
pass(1) = ( norm(A - X) < 100*tol); 

% check the @chebfun2/chebpolyval2 command.
f = chebfun2(@(x,y) cos(x.*y) ); 
[A1, A2, A3] = chebpolyval2( f ); 
X = chebpolyval2( f ); 
pass(2) = ( norm(X - A1*A2*A3.') < tol);

% Check when degree in x \neq degree in y: 
A = rand(10,11);
f = chebfun2(A);
B = chebpolyval2(f);
pass(2) = ( norm(A - B) < tol);
end