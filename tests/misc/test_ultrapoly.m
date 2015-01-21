% Test file for ULTRAPOLY

function pass = test_ultrapoly(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

tol = 1e-14;

lambda = 2.1; 
U = ultrapoly(10, lambda);
r = roots( U ); 
exact = jacpts( 10, lambda-.5, lambda -.5);
pass(1) = norm( sort(r) - sort(exact) ) < tol ; 

lambda = 1.9; 
U = ultrapoly(11, lambda);
r = roots( U ); 
exact = jacpts( 11, lambda-.5, lambda -.5);
pass(2) = norm( sort(r) - sort(exact) ) < tol ; 

lambda = 1.9; 
U = ultrapoly(0:15, lambda);
J = jacpoly(0:15, lambda-.5, lambda-.5);
x = linspace(-1,1,100)';
UU = U(x,:); 
JJ = J(x,:); 
A = UU./JJ; 
pass(3) = norm( A - ones(100,1)*A(1,:) ) < 1e4*tol;

lambda = 4; 
U = ultrapoly(0:15, lambda);
J = jacpoly(0:15, lambda-.5, lambda-.5);
x = linspace(-1,1,100)';
UU = U(x,:); 
JJ = J(x,:); 
A = UU./JJ; 
pass(4) = norm( A - ones(100,1)*A(1,:) ) < 1e5*tol;

lambda = .25;
U = ultrapoly([0 2 3], lambda);
C = chebcoeffs(U);
tru = [1 0.0625 0
       0 0      0.078125
  	   0 0.3125 0
       0 0      0.234375];
pass(5) = ( all(size(C) == size(tru)) ) && ( norm(C(:) - tru(:)) < tol );

end