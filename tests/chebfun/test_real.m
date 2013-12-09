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
pass(1) = isempty(real(f));

f = chebfun(@(x) exp(1i*x).*abs(x - 0.1), [-1 1], pref);
ref = real(f);
ref_exact = @(x) cos(x).*abs(x - 0.1);
err = feval(ref, xr) - ref_exact(xr);
pass(2) = isreal(ref) && (norm(err(:), inf) < 10*ref.vscale.*ref.epslevel);

f = chebfun(@(x) [exp(1i*x).*abs(x - 0.1) 1i*exp(-1i*x).*cos(x)], [-1 1], pref);
ref = real(f);
ref_exact = @(x) [cos(x).*abs(x - 0.1) sin(x).*cos(x)];
err = feval(ref, xr) - ref_exact(xr);
pass(3) = isreal(ref) && (norm(err(:), inf) < 10*ref.vscale.*ref.epslevel);

f = chebfun(@(x) 1i*cos(x).*abs(x - 0.1), [-1 1], pref);
ref = real(f);
pass(4) = iszero(ref);

% Integration of SINGFUN:
dom = [-2 7];
pow = -1.64;
f = chebfun(@(x) (x-dom(1)).^pow, dom, 'exps', [pow 0], 'splitting', 'on');
pass(5) = ( isequal(real(f), f) );

end
