function pass = test_chebcoeffs3( pref ) 
% Test chebcoeffs3

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb3Prefs.chebfun3eps;

% Create a rank-(3,3,3) function:
m = 8;
n = 10;
p = 5;
Tm = chebpoly(m);
Tn = chebpoly(n);
Tp = chebpoly(p);
f = outerProd(Tm, Tm, Tm) + outerProd(Tn, Tn, Tn) + outerProd(Tp, Tp, Tp);

% Create expected tensor of coefficients for f:
fCoeffs = chebcoeffs3(f);
coeffsExact = zeros(m+1, n+1, p+1); 
coeffsExact(m+1, m+1, m+1) = 1; 
coeffsExact(n+1, n+1, n+1) = 1; 
coeffsExact(p+1, p+1, p+1) = 1; 
pass(1) = norm( fCoeffs(:) - coeffsExact(:) ) < tol;

% check the reverse:
C = chebfun3.vals2coeffs( chebpolyval3(f) );
pass(2) = norm(C(:) - coeffsExact(:)) < tol;

end