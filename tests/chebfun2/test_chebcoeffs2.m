function pass = test_chebcoeffs2( pref ) 
% Test chebpoly2 

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb2Prefs.chebfun2eps;
j = 1; 

% Rank-2 function
n = 10;
m = 8;
Tn = chebpoly(n);
Tm = chebpoly(m);
f = Tn * Tn' + Tm * Tn'; 
X = chebcoeffs2( f );
Exact = zeros(n+1); Exact(n+1,n+1) = 1; Exact(m+1,n+1) = 1; 
pass(j) = norm( X - Exact ) < tol; j = j + 1; 

%f = Tm * Tn';
% check inverses
Z = chebfun2.vals2coeffs( chebpolyval2( f ) );
pass(j) = norm( Z - Exact ) < tol; j = j + 1; 

end