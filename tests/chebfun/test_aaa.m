function pass = test_aaa(pref)
% Test for aaa.  This has been written to be independent of Chebfun.  

% Get preferences.
if ( nargin < 1 )
    tol = 1e4*eps;
end

warning('off', 'AAA:Froissart');

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

% Two very short cases.
Z = [0, 1];
F = [1, 2];
r = aaa(F, Z);
pass(10) = ( norm(F - r(Z), inf) < tol );
Z = [0, 1, 2];
F = [1, 0, 0];
r = aaa(F, Z);
pass(11) = ( norm(F - r(Z), inf) < tol );

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
Z = randn(10000,1)+3i*randn(10000,1);
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

% Make sure Lawson matches minimax and degree differs from mmax
f = @(x) exp(x);
xx = linspace(-1,1);
err_minimax = 1.550669058714149e-7;
r = aaa(f,'degree',3); err_aaa = norm(f(xx) - r(xx),inf);
pass(21) = (err_aaa/err_minimax < 1.1);
r = aaa(f,'mmax',4); err_aaa = norm(f(xx) - r(xx),inf);
pass(22) = (err_aaa/err_minimax > 1.1);

% Make sure Lawson bails out if unsuccessful because of machine precision
xx = linspace(-1,1);
r = aaa(@tanh,xx); err1 = norm(tanh(xx) - r(xx),inf);
r = aaa(@tanh,xx,'mmax',40); err2 = norm(tanh(xx) - r(xx),inf);
pass(23) = abs(err2/err1 - 1) < 1.01; 

% Make sure Lawson bails out if unsuccessful because of symmetry
Z = exp(2i*pi*(1:500)'/500); F = log(2-Z.^4); n = 15;
r = aaa(F,Z,'mmax',n+1,'lawson',0); err1 = norm(F - r(Z),inf);
r = aaa(F,Z,'mmax',n+1); err2 = norm(F - r(Z),inf);
pass(24) = abs(err2/err1 - 1) < 1.01; 

% Make sure Lawson bails out if unsuccessful because of troublesome poles
Za = linspace(-3,-1,1000)'; Zb = linspace(1,3,1000)'; Z = [Za; Zb];
F = [sign(Za); sign(Zb)]; n = 12; 
r = aaa(F,Z,'mmax',n+1,'lawson',0); err1 = norm(F - r(Z),inf);
r = aaa(F,Z,'mmax',n+1); err2 = norm(F - r(Z),inf);
pass(25) = abs(err2/err1 - 1) < 1.01; 

% Degree option 
Z = linspace(-1, 1, 1000);
F = exp(Z);
[r, pol] = aaa(F, Z, 'degree', 3);
[r2, pol2] = aaa(F, Z, 'mmax', 4);
pass(26) = (numel(pol) == 3);
pass(27) = (numel(pol2) == 3);

% Bug reported by Williams Johns in issue 2423
X = [1 2 3];
F = [1 0 0];
r = aaa(F,X);
pass(28) = (norm(F-r(X)) == 0);

X = chebpts(200); F = max(X,0);
deg = 4;
r = aaa(F,X,'degree',4,'lawson',100,'damping',0.2);
err = norm(F-r(X),inf); pass(29) = abs(err-.006) < .002;
r = aaa(F,X,'degree',4,'lawson',100,'damping',0.5,'sign',1);
err = norm(F-r(X),inf); pass(30) = abs(err-.006) < .002;

X = chebpts(1000,[0 10]); F = 1./(1+exp(5./(X-2)));   % Fermi-Dirac
r = aaa(F,X,'degree',12,'lawson',100,'damping',0.85,'sign',1);
err = norm(F-r(X),inf); pass(31) = abs(err-.000035) < .0001;

f = @(x) max(x,0);
xx = linspace(-1,1,300);
r = aaa(f(xx),xx,'degree',8,'damping',.5,'lawson',200);
% r = aaa(f,'degree',8,'damping',.5,'lawson',200);
% xx = linspace(-1,1,300);
err = norm(f(xx)-r(xx),inf); pass(32) = abs(err-.0006) < .001;

Z = linspace(-1,1,100); F = exp(Z);
r = aaa(F, Z, 'deriv_deg', 1);
err = norm(F-r{2}(Z),inf); pass(33) = abs(err) < 1e-12;

Z = linspace(-1,1,100); F = sin(Z);
r = aaa(F, Z, 'deriv_deg', 4);
pass(34) = iscell(r) && isequal(size(r), [1, 5]);

Z = linspace(-1,1,100); F = sin(Z);
r = aaa(F, Z, 'deriv_deg', 4);
err = norm(F-r{5}(Z),inf);
pass(35) = abs(err) < 1e-7;

Z = linspace(-1,1,100); F = Z.^3;
r = aaa(F, Z, 'deriv_deg', 1);
err = norm(3*Z.^2-r{2}(Z),inf);
pass(36) = abs(err) < 1e-12;

Z = linspace(-1,1,100); F = 1e-100*exp(Z);
r = aaa(F, Z, 'deriv_deg', 1);
err = norm(F-r{2}(Z),inf); pass(37) = abs(err) < 1e-112;

Z = linspace(-1e100,1e100,100); F = exp(1e-100*Z);
r = aaa(F, Z, 'deriv_deg', 1);
err = norm(1e-100*F-r{2}(Z),inf); pass(38) = abs(err) < 1e-112;

Z = exp(2i*pi*(0:99)/100); F = sqrt(2-Z);
rr = aaa(F, Z, 'degree', 6, 'deriv_deg', 1);
r = rr{1};
err = abs( norm(F-r(Z),inf) - 1.19e-11 ); pass(39) = abs(err) < 1e-12;
rp = rr{2};
pass(40) = abs(rp(1)+0.5) < 1e-6;

Z = linspace(-1,1,100); F = Z + 1./(Z-1.5);
[r,pol,res] = aaa(F, Z); % residual computation check
[~,ii] = min(abs(pol-1.5));
err = res(ii)-1;
pass(41) = abs(err) < 1e-8;

Z = logspace(-15,0,300)'; F = Z.^(1/2); 
r = aaa(F, Z); % this test checks diag scaling is working
ZZ = logspace(-15,0,500)';
err = norm(r(ZZ)-sqrt(ZZ),inf);
pass(42) = err < 1e-8;

%%% Noise chop tests

% Test basic behavior: chopping occurs and consistent outputs 
state = rng;
rng(0)
X = linspace(-1,1,500);
F = sin(10*X) + 1e-8*randn(1,500);
[~,~,~,~,zj0,fj0,~,errvec_full,~,svals0] = aaa(F, X);
[~,pol,res,zer,zj,fj,wj,errvec_chop,~,svals] = ...
    aaa(F, X, 'noise_chop', 1);
pass(43) = isvector(errvec_chop) && isvector(errvec_full);
pass(44) = length(errvec_full) == length(svals0);
pass(45) = length(errvec_chop) < length(errvec_full);
n = length(errvec_chop);
pass(46) = all(cellfun(@length, ...
    {zj, fj, wj, svals}) == n);
m = length(pol);
[~,poln,resn,zern,zjn,fjn,wjn,errvecn,~,svalsn] = aaa(F, X, 'mmax', n);
pass(47) = (m == length(res)) && (length(zer) <= n-1) && ...
    isequal(zj0(1:n), zj) && isequal(fj0(1:n), fj) && ...
    isequal(errvec_full(1:n), errvec_chop) && ...
    isequal(svals0(1:n), svals) && isequal(zjn, zj) && ...
    isequal(fjn, fj) && isequal(wjn, wj) && ...
    isequal(errvecn, errvec_chop) && isequal(svalsn, svals) && ...
    isequal(poln, pol) && isequal(resn, res) && isequal(zern, zer);
[~,~,~,~,~,~,~,errvec_off] = aaa(F, X, 'noise_chop', 0);
pass(48) = isequal(errvec_full, errvec_off);

% Test that noise chop does not affect clean data
Fclean = sin(10*X);
[~,~,~,~,~,~,~,errvec_clean] = aaa(Fclean, X);
[~,~,~,~,~,~,~,errvec_clean_full] = aaa(Fclean, X, 'noise_chop', 0);
pass(49) = isequal(errvec_clean, errvec_clean_full);

% Test logical values for options 
[~,~,~,~,~,~,~,errvec_true] = aaa(Fclean, X, 'noise_chop', true);
[~,~,~,~,~,~,~,errvec_false] = aaa(Fclean, X, 'noise_chop', false);
[~,~,~,~,~,~,~,errvec_zero] = aaa(Fclean, X, 'noise_chop', 0);
[~,~,~,~,~,~,~,errvec_one] = aaa(Fclean, X, 'noise_chop', 1);
pass(50) = isequal(errvec_true, errvec_one);
pass(51) = isequal(errvec_false, errvec_zero);
pass(52) = isequal(errvec_false, errvec_clean);

% Test on short data with noise chop on and off
rng(0)
Fshort = sin(10*X) + 1e-8*randn(1,500);
[~,~,~,~,~,~,~,errvec_short] = aaa(Fshort, X, 'mmax', 20);
rng(0)
Fshort = sin(10*X) + 1e-8*randn(1,500);
[~,~,~,~,~,~,~,errvec_short_full] = aaa(Fshort, X, 'mmax', 20, ...
    'noise_chop', 0);
pass(53) = isequal(errvec_short, errvec_short_full);

Xlawson = linspace(-1,1,200);
Flawson = max(Xlawson,0);
[r0,~,~,~,zj0,~,~,errvec0] = aaa(Flawson, Xlawson, 'degree', 4, ...
    'lawson', 100, 'damping', 0.2);
[r1,~,~,~,zj1,~,~,errvec1] = aaa(Flawson, Xlawson, 'degree', 4, ...
    'lawson', 100, 'damping', 0.2, 'noise_chop', 1);
pass(54) = isequal(errvec0, errvec1) && length(zj0) == length(zj1) && ...
    norm(r0(Xlawson)-r1(Xlawson), inf) == 0;

[r0,~,~,~,zj0,~,~,errvec0] = aaa(Flawson, Xlawson, 'degree', 4, ...
    'lawson', 100, 'damping', 0.5, 'sign', 1);
[r1,~,~,~,zj1,~,~,errvec1] = aaa(Flawson, Xlawson, 'degree', 4, ...
    'lawson', 100, 'damping', 0.5, 'sign', 1, 'noise_chop', 1);
pass(55) = isequal(errvec0, errvec1) && length(zj0) == length(zj1) && ...
    norm(r0(Xlawson)-r1(Xlawson), inf) == 0;

[r_auto,~,~,~,~,~,~,errvec_auto,wt_auto,svals_auto] = aaa(@exp, ...
    'noise_chop', 0);
[r_auto_chop,~,~,~,~,~,~,errvec_auto_chop,~,svals_auto_chop] = ...
    aaa(@exp, 'noise_chop', 1);
xx = linspace(-1,1);
pass(56) = isa(r_auto, 'function_handle');
pass(57) = ~isempty(errvec_auto) && (numel(errvec_auto) == numel(svals_auto));
pass(58) = isvector(wt_auto) && all(isnan(wt_auto));
pass(59) = isequal(errvec_auto, errvec_auto_chop) && ...
    isequal(svals_auto, svals_auto_chop) && ...
    norm(r_auto(xx)-r_auto_chop(xx), inf) == 0;

rng(0)
F = tanh(20*X) + 1e-8*randn(1,500);
[r_lawson_chop,~,~,~,zj_lawson_chop,~,~,errvec_lawson_chop] = ...
    aaa(F, X, 'noise_chop', 1, 'lawson', 5, 'sign', 1);
rng(0)
F = tanh(20*X) + 1e-8*randn(1,500);
[~,~,~,~,zj_lawson_full,~,~,errvec_lawson_full] = ...
    aaa(F, X, 'noise_chop', 0, 'lawson', 5, 'sign', 1);
pass(60) = length(errvec_lawson_chop) < length(errvec_lawson_full);
pass(61) = length(zj_lawson_chop) < length(zj_lawson_full);
pass(62) = all(isfinite(r_lawson_chop(X))) && ...
    norm(F-r_lawson_chop(X), inf) < 1e-5;

rng(0)
X = linspace(-1,1,80);
F = sin(10*X) + 1e-8*randn(size(X));
[~,~,~,~,~,~,~,e0] = aaa(F,X,'mmax',80,'noise_chop',0);
[~,~,~,~,~,~,~,e1,~,s1] = aaa(F,X,'mmax',80,'noise_chop',1);
pass(63) = length(e1) == length(s1);
pass(64) = length(e1) <= floor(length(X)/2);
pass(65) = length(e0) > length(e1);

rng(state)


warning('on', 'AAA:Froissart');

end
