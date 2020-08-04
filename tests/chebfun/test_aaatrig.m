function pass = test_aaatrig(pref)
% Test for aaatrig().

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end
tol = 1e4*pref.chebfuneps;

warning('off', 'CHEBFUN:aaatrig:Froissart');

% Do tests for both odd and even approximants
for j = 1:2

if j == 1
    form = 'odd';
else 
    form = 'even';
end

Z = linspace(0, 2*pi, 1000);
F = exp(cos(Z));
[r, ~, ~, ~, zj] = aaatrig(F, Z, form);
pass(1,j) = ( norm(F - r(Z), inf) < tol );
pass(2,j) = isnan(r(nan));                  % check that r(NaN) = NaN
pass(3,j) = isnan(r(inf));                  % real infinite evaluation should return nan
pass(4,j) = all(~isinf(r([1i;-1i]*inf)));   % imaginary infinite evaluation should returns a finite number
m1 = length(zj);
[r, ~, ~, ~, zj] = aaatrig(F, Z, form, 'mmax', m1 - 1);
pass(4,j) = ( length(zj) == m1 - 1 );
[r, ~, ~, ~, zj] = aaatrig(F, Z, form, 'tol', 1e-3);
pass(5,j) = ( length(zj) < m1 );

Z = linspace(0, 2*pi, 1000);
F = @(z) sin(3*z).*cos(7*z)./sin(2*z-1i);
[r, pol, res, zer] = aaatrig(F, Z, form,'cleanup',0);
pass(6,j) = ( norm(F(Z) - r(Z), inf) < 10*tol );
pass(7,j) = ( min( abs(zer-pi/2)) < tol );
pass(8,j) = (min( abs(real(pol - 0.5i - pi/2))) < tol );

% Case |Z| = 2: needs special treatment.
Z = [0, 1];
F = [1, 2];
r = aaatrig(F, Z, form);
pass(9,j) = ( norm(F - r(Z), inf) < tol );

% Check for exact scale-invariance
Z = linspace(0.3,1.5);
F = exp(Z)./(1+1i);
r1 = aaatrig(F, Z, form);
r2 = aaatrig(2^311*F, Z, form);
r3 = aaatrig(2^-311*F, Z, form);
pass(10,j) = ( r1(0.2i) == 2^-311*r2(0.2i) ); 
pass(11,j) = ( r1(1.4) == 2^311*r3(1.4) ); 

rng(0); Z = pi*randn(10000,1)+1i*randn(10000,1);
f = @(z) log(5-tan(z/2))./(1+tan(z/2).^2);
r = aaatrig(f(Z),Z,form);
pass(12,j) = ( abs(r(0) - f(0)) < tol );
 
% Test behavior for string inputs:
Z = linspace(-1,1,10001);
r1 = aaatrig(@(x) abs(tan(x/2)), Z, form);
r2 = aaatrig('abs(tan(x/2))', Z, form);
x = -1 +2*rand(1);
pass(13,j) = ( r1(x) == r2(x) );

% Test that constructor does not fail when a data value is infinite:
Z = linspace(-1,1);
r = aaatrig(gamma(Z), Z, form, 'cleanup',0);
pass(14,j) = ( abs(r(0.63) - gamma(0.63)) < 1e-3 );

% Test for NaNs
X = linspace(0,pi);
F = sin(X)./X;
r = aaatrig(F,X,form);
pass(15,j) = ( abs(r(2) - sin(2)/2) < 1e-3 );
 
% A couple of tests of residues
X = linspace(-1.337,2,537);
[r,pol,res] = aaatrig(cot(X/2).*exp(1i*tan(X/2)), X);
ii = find(abs(mod(real(pol),2*pi))<1e-8);
pass(16,j) = abs(res(ii)-2) < 1e-10;
[r,pol,res] = aaatrig(tan(X/2).*log(5+sin(X)),X);
ii = find(abs(pol - pi)<1e-8);
pass(17,j) = abs(res(ii) + 2*log(5)) < 1e-8;

% Make sure Lawson bails out if unsuccessful because of machine precision
xx = linspace(-1,1);
r = aaatrig(@tanh,xx,form); err1 = norm(tanh(xx) - r(xx),inf);
r = aaatrig(@tanh,xx,form,'mmax',40); err2 = norm(tanh(xx) - r(xx),inf);
pass(18,j) = abs(err2/err1 - 1) < 1.01; 
 
% Make sure Lawson bails out if unsuccessful because of symmetry
Z = exp(2i*pi*(1:500)'/500); F = log(2-sin(Z).^4); n = 15;
r = aaatrig(F,Z,form,'mmax',n+1,'lawson',0); err1 = norm(F - r(Z),inf);
r = aaatrig(F,Z,form,'mmax',n+1); err2 = norm(F - r(Z),inf);
pass(19,j) = abs(err2/err1 - 1) < 1.01; 

% Make sure Lawson bails out if unsuccessful because of troublesome poles
Za = chebpts(1000,[-3 -1]); Zb = chebpts(1000,[1,3]); Z = [Za; Zb];
F = [sign(Za); sign(Zb)]; n = 12; 
r = aaatrig(F,Z, form, 'mmax',n+1,'lawson',0); err1 = norm(F - r(Z),inf);
r = aaatrig(F,Z, form, 'mmax',n+1); err2 = norm(F - r(Z),inf);
pass(20,j) = abs(err2/err1 - 1) < 1.01; 

% Degree option 
Z = linspace(-1, 1, 1000);
F = exp(Z);
[r, pol] = aaatrig(F, Z, form, 'degree', 3);
[r2, pol2] = aaatrig(F, Z, form, 'mmax', 4);
pass(21,j) = (numel(pol) <= 4);
pass(22,j) = (numel(pol2) <= 4);
pass(23,j) = (norm(r(Z)-r2(Z)) < 1e-10);

end

warning('on', 'CHEBFUN:aaatrig:Froissart');

pass = pass(:).';

end
