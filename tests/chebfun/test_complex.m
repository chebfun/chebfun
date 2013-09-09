% Test file for @chebfun/complex.m.

function pass = test_complex(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfun.pref;
end

pref.chebfun.splitting = 1;

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Check empty cases.
f = chebfun(@sin, [-1 1], pref);
pass(1) = isempty(complex(chebfun(), chebfun()));
pass(2) = isempty(complex(f, chebfun()));
pass(3) = isempty(complex(chebfun(), f()));

% Check a known result.
ref = chebfun(@cos, [-1 0 1], pref);
imf = chebfun(@sin, [-1 0 1], pref);
f = complex(ref, imf);
f_exact = @(x) exp(1i*x);
err = feval(f, xr) - f_exact(xr);
pass(4) = norm(err, inf) < 10*f.vscale*f.epslevel;

% Check behavior for an array-valued function.
ref = chebfun(@(x) [cos(x) -sin(x)], [-1 0 1], pref);
imf = chebfun(@(x) [sin(x) cos(x)], [-1 0 1], pref);
f = complex(ref, imf);
f_exact = @(x) [exp(1i*x) 1i*exp(1i*x)];
err = feval(f, xr) - f_exact(xr);
pass(5) = norm(err(:), inf) < 10*f.vscale*f.epslevel;

% Check forming from just a single real chebfun.
pass(6) = isequal(ref, complex(ref));

% Check forming from a chebfun and a real scalar.
alpha = -0.194758928283640;

f = complex(ref, alpha);
f_exact = @(x) [cos(x) -sin(x)] + alpha*1i;
err = feval(f, xr) - f_exact(xr);
pass(7) = norm(err(:), inf) < 10*f.vscale*f.epslevel;

f = complex(alpha, imf);
f_exact = @(x) alpha + 1i*[sin(x) cos(x)];
err = feval(f, xr) - f_exact(xr);
pass(8) = norm(err(:), inf) < 10*f.vscale*f.epslevel;

% Check error conditions.
f = chebfun(@sin, [-1 1]);
g = chebfun(@(x) exp(1i*x), [-1 1]);

try
    h = complex(g, f);
    pass(9) = false;
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:complex:AisNotReal');
end

try
    h = complex(f, g);
    pass(10) = false;
catch ME
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:complex:BisNotReal');
end

end
