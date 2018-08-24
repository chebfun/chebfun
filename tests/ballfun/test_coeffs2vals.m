function pass = test_coeffs2vals( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps; 

% Test with function sin(r)*th*lam
f = ballfun(@(x,y,z)sin(x).*y.*z, 'cart' );

exact = f.coeffs; 
cfs = ballfun.vals2coeffs(ballfun.coeffs2vals(f.coeffs));

pass(1) = norm( exact(:) - cfs(:), inf) < tol; 

end
