function pass = test_aaa(pref)
% Test for aaa().

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end
tol = 1e+4 * pref.chebfuneps;


Z = linspace(-1, 1, 1000);
F = exp(Z);
[r, ~, ~, ~, zj] = aaa(F, Z);
pass(1) = ( norm(F - r(Z), inf) < tol );
pass(2) = isnan(r(nan));                        % check that r(NaN) = NaN
pass(3) = ~isinf(r(inf));                       % r(inf) = sum(w.*f)/sum(w)
m1 = length(zj);
[r, ~, ~, ~, zj] = aaa(F, Z, 'mmax', m1 - 1);
pass(4) = ( length(zj) == m1 - 1 );
[r, ~, ~, ~, zj] = aaa(F, Z, 'tol', 1e-3);
pass(5) = ( length(zj) < m1 );


%
Z = linspace(-1, 1, 1000);
F = @(z) tan(pi*z);
[r, pol, res, zer] = aaa(F, Z);
pass(6) = ( norm(F(Z) - r(Z), inf) < 10*tol );
pass(7) = ( min(abs(zer)) < tol );
pass(8) = ( min(abs(pol - 0.5)) < tol );
pass(9) = ( min(abs(res)) > 1e-13 );        % Test for spurious poles.


% Case |Z| = 2: needs special treatment.
Z = [0, 1];
F = [1, 2];
r = aaa(F, Z);
pass(10) = ( norm(F - r(Z), inf) < tol );
pass(11) = ( r(inf) == -inf );

end
