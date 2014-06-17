function pass = test_divgrad( pref ) 
% Test DIVGRAD 
if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 50*pref.eps; 

% Check definition: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y)); 
divgradF = diffx(F(1), 2) + diffy(F(2), 2); 
pass(1) = ( norm(divgradF - divgrad(F)) < tol );

end