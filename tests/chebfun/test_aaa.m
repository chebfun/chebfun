function pass = test_aaa(pref)
% Test for aaa().  Many new tests added August 2019.

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end
tol = 1e4*pref.chebfuneps;

warning('off', 'CHEBFUN:aaa:Froissart');

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

% Test that constructor does not fail when a data value is infinite:
Z = linspace(-1,1);
r = aaa(gamma(Z),Z);
pass(17) = ( abs(r(0.63) - gamma(0.63)) < 1e-3 );

% Test for NaNs
X = linspace(0,20);
F = sin(X)./X;
r = aaa(F,X);
pass(18) = ( abs(r(2) - sin(2)/2) < 1e-3 );

% A couple of tests of residues
X = linspace(-1.337,2,537);
[r,pol,res] = aaa(exp(X)./X, X);
ii = find(abs(pol)<1e-8);
pass(19) = abs(res(ii)-1) < 1e-10;
[r,pol,res] = aaa((1+1i)*gamma(X),X);
ii = find(abs(pol-(-1))<1e-8);
pass(20) = abs(res(ii)+(1+1i)) < 1e-10;

% Make sure Lawson matches minimax:
x = chebfun('x'); f = exp(x);
xx = linspace(-1,1);
[p,q,r] = minimax(f,3,3); err_minimax = norm(f(xx) - r(xx),inf);
r = aaa(f,'mmax',4); err_aaa = norm(f(xx) - r(xx),inf);
pass(21) = (err_aaa/err_minimax < 1.1);

% Make sure Lawson bails out if unsuccessful because of machine precision
xx = linspace(-1,1);
r = aaa(@tanh,xx); err1 = norm(tanh(xx) - r(xx),inf);
r = aaa(@tanh,xx,'mmax',40); err2 = norm(tanh(xx) - r(xx),inf);
pass(22) = abs(err2/err1 - 1) < 1.01; 

% Make sure Lawson bails out if unsuccessful because of symmetry
Z = exp(2i*pi*(1:500)'/500); F = log(2-Z.^4); n = 15;
r = aaa(F,Z,'mmax',n+1,'lawson',0); err1 = norm(F - r(Z),inf);
r = aaa(F,Z,'mmax',n+1); err2 = norm(F - r(Z),inf);
pass(23) = abs(err2/err1 - 1) < 1.01; 

% Make sure Lawson bails out if unsuccessful because of troublesome poles
Za = chebpts(1000,[-3 -1]); Zb = chebpts(1000,[1,3]); Z = [Za; Zb];
F = [sign(Za); sign(Zb)]; n = 12; 
r = aaa(F,Z,'mmax',n+1,'lawson',0); err1 = norm(F - r(Z),inf);
r = aaa(F,Z,'mmax',n+1); err2 = norm(F - r(Z),inf);
pass(24) = abs(err2/err1 - 1) < 1.01; 

% Degree option 
Z = linspace(-1, 1, 1000);
F = exp(Z);
[r, pol] = aaa(F, Z, 'degree', 3);
[r2, pol2] = aaa(F, Z, 'mmax', 4);
pass(25) = (numel(pol) == 3);
pass(26) = (numel(pol2) == 3);
pass(27) = (norm(r(Z)-r2(Z)) < 1e-10);

warning('on', 'CHEBFUN:aaa:Froissart');

end
