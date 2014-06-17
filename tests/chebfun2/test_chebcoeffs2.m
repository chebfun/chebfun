function pass = test_chebpoly2( pref ) 
% Test chebpoly2 

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.eps; 
j = 1; 

% Rank-2 function
T10 = chebpoly(10);
T8 = chebpoly(8);
f = T10 * T10' + T8 * T10'; 
X = chebcoeffs2( f );
Exact = zeros(11); Exact(1,1) = 1; Exact(3,1) = 1; 
pass(j) = norm( X - Exact ) < tol; j = j + 1; 

% check inverses
Z = chebfun2.vals2coeffs( chebpolyval2( f ) );
pass(j) = norm( Z - Exact ) < tol; j = j + 1; 

end