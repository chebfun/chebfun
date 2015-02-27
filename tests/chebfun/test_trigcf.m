% Test for trigremez.m.
function pass = test_trigcf(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Set an error tolerance.
tol = 1.0e-12;

%% check for emptiness:
p = trigcf(chebfun);
pass(1) = isempty(p);


%% check for a simple trigfun 
fh = @(x) cos(4*pi*x) + 1;
f = chebfun(fh, [0, 2], 'trig');
[p, s] = trigcf(f, 4);
pass(2) = norm(p-f, inf) < tol && abs(s) < tol;

%%
fh = @(x) exp(sin(pi*x));
f = chebfun(fh, [0, 2], 'trig');
n = 12;
pcf = trigcf(f, n);
[prm, err, status] = trigremez(f, n);
%%
errCF = chebfun(pcf-f);
errRM = chebfun(prm-f);
semilogy(abs(errCF))
hold on
semilogy(abs(errRM), 'r')
hold off

%%

%%
err = f - p;
pass(1) = isPeriodicTech(p) & (length(p) == 2*n+1);
pass(2) = norm(abs(err(status.xk)) - errMax, inf) < 100*tol;
pass(3) = norm(p-1,inf) < 100*tol;
n = 4;

%% Pass a non-trig chebfun and reproduce it:
f = chebfun(f);
[p, errMax, status] = trigremez(f, n);



%%
% Test trigonometric polynomail cases
f = chebfun(@(x) exp(cos(pi*x) + sin(2*pi*x)), 'trig');
M = ceil((length(f)-1)/2);
pass(9) = norm(cf(f,M) - f, inf) < 100*eps*f.vscale;
pass(10) = norm(cf(f, M+10) - f, inf) < 100*eps*f.vscale;
f = chebfun(@(x) cos(4*pi*x), 'trig');
[p, ~, ~, s] = cf(f, 3);
pass(11) = norm(p, inf) < 100*eps*f.vscale;
pass(12) = abs(s - 1) < 100*eps*f.vscale;