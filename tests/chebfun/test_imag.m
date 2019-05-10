% Test file for @chebfun/real.m.

function pass = test_imag(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

pref.splitting = 1;

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Run a few straightforward tests.
f = chebfun();
pass(1) = isempty(imag(f));

f = chebfun(@(x) exp(1i*x).*abs(x - 0.1), [-1 1], pref);
imf = imag(f);
imf_exact = @(x) sin(x).*abs(x - 0.1);
err = feval(imf, xr) - imf_exact(xr);
pass(2) = isreal(imf) && (norm(err(:), inf) < 10*vscale(imf)*eps);

f = chebfun(@(x) [exp(1i*x).*abs(x - 0.1) 1i*exp(-1i*x).*cos(x)], [-1 1], pref);
imf = imag(f);
imf_exact = @(x) [sin(x).*abs(x - 0.1) cos(x).^2];
err = feval(imf, xr) - imf_exact(xr);
pass(3) = isreal(imf) && (norm(err(:), inf) < 1e2*vscale(imf)*eps);
    

f = chebfun(@(x) 1i*cos(x).*abs(x - 0.1), [-1 1], pref);
imf = real(f);
pass(4) = iszero(imf);

%% Test on singular function:

dom = [-2 7];
pow = -1.64;
f = chebfun(@(x) 1i*(x-dom(1)).^pow, dom, 'exps', [pow 0], 'splitting', 'on');
pass(5) = ( isequal(imag(f), -1i*f) );

%% Test on function defined on unbounded domain:

% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];
domCheck = [-1e6 -3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% Array-valued function:
op = @(x) x.*exp(x) + 1i*(1-exp(x))./x;
opg = @(x) (1-exp(x))./x;
f = chebfun(op, dom);
g = imag(f);
gVals = feval(g, x);
gExact = opg(x);
err = gVals - gExact;
pass(6) = norm(err, inf) < 3e1*eps*vscale(f);

end
