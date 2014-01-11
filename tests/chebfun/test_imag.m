% Test file for @chebfun/real.m.

function pass = test_real(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

pref.enableBreakpointDetection = 1;

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
pass(2) = isreal(imf) && (norm(err(:), inf) < 10*imf.vscale.*imf.epslevel);

f = chebfun(@(x) [exp(1i*x).*abs(x - 0.1) 1i*exp(-1i*x).*cos(x)], [-1 1], pref);
imf = imag(f);
imf_exact = @(x) [sin(x).*abs(x - 0.1) cos(x).^2];
err = feval(imf, xr) - imf_exact(xr);
pass(3) = isreal(imf) && (norm(err(:), inf) < 10*imf.vscale.*imf.epslevel);

f = chebfun(@(x) 1i*cos(x).*abs(x - 0.1), [-1 1], pref);
imf = real(f);
pass(4) = iszero(imf);

end
