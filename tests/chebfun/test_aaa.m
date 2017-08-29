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

% Check for exact scale-invariance
Z = linspace(0.3,1.5);
F = exp(Z)./(1+1i);
r1 = aaa(F, Z);
r2 = aaa(2^311*F, Z);
r3 = aaa(2^-311*F, Z);
pass(12) = ( r1(0.2i) == 2^-311*r2(0.2i) ); 
pass(13) = ( r1(1.4) == 2^311*r3(1.4) ); 

% Make sure the gamma function gives something reasonable:
r = aaa(@gamma);
pass(14) = ( abs(r(1.5) - gamma(1.5)) < 1e-3 );

%
rng(0); Z = randn(10000,1)+3i*randn(10000,1);
f = @(z) log(5-z)./(1+z.^2);
r = aaa(f(Z),Z);
pass(15) = ( abs(r(0) - f(0)) < tol );

% Test behavior for string inputs:
Z = linspace(-1,1,10001);
r1 = aaa(@(x) abs(x), Z);
r2 = aaa('abs(x)', Z);
x = -1 +2*rand(1);
pass(16) = ( r1(x) == r2(x) );

end
