function pass = test_coeffs2vals( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps; 

% Test with function sin(r)*th*lam
f = ballfun(@(x,y,z)sin(x).*y.*z);

exact = f.coeffs; 
cfs = ballfun.vals2coeffs(ballfun.coeffs2vals(f.coeffs));

pass(1) = norm( exact(:) - cfs(:), inf) < tol; 

% Example 2
f = @(lam)exp(1i*lam);
lam = pi*trigpts(3);
vals = f(lam).';
cfs = [0,0,1];
pass(2) = norm(ballfun.coeffs2vals(cfs)-vals) < tol;

if (nargout > 0)
    pass = all(pass(:));
end

end
