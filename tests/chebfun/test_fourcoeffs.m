function pass = test_fourcoeffs(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

% Test on simple combination of Fourier modes on the standard
% domain of [-pi pi].
f_test = @(x) 2 + 2*cos(2*x) + sin(x);
dom = [-pi pi];
f = chebfun(f_test,dom,'periodic');
c = fourcoeffs(f);
c_exact = [1 -0.5i 2 0.5i 1];
err = c-c_exact;
pass(1) = norm(err,inf) < vscale(f).*epslevel(f);

% Now change domains and check that the same result is given for the 
% coefficients.
dom = [0.1 0.1+2*pi];
f = chebfun(f_test,dom,'periodic');
c = fourcoeffs(f);
err = c-c_exact;
pass(2) = norm(err,inf) < vscale(f).*epslevel(f);

end