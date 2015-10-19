% Test file for chebfun/circconv

function pass = test_circconv(pref)

% Grab some preferences.
if ( nargin == 0 )
    pref = chebfunpref();
end

testclass = trigtech();

% Generate a random point to use as test values.
seedRNG(6178);
x = 2 * rand - 1;

% Check operation in the face of empty arguments.
f = chebfun();
g = chebfun();
pass(1) = (isempty(circconv(f, g)) && isempty(circconv(g, f)));

% No support for quasimatrices.
f = chebfun(@(x) [cos(10*x) sin(6*x)], [-pi pi], 'periodic'); 
try
    h = circconv(f, f);
catch ME
    pass(2) = strcmpi(ME.message, 'No support for array-valued CHEBFUN objects.');
end

% No support for Chebyshev-based CHEBFUN objects.
f = chebfun(@(x) sin(20*x), [-pi pi]); 
g = chebfun(@(x) cos(cos(x)), [-pi pi], 'periodic');
try
    h = circconv(f, g);
catch ME
    pass(3) = strcmpi(ME.message, 'No support for Chebyshev-based CHEBFUN objects.');
end

% Check domains.
f = chebfun(@(x) sin(20*x), [-2*pi pi], 'periodic'); 
g = chebfun(@(x) cos(cos(x)), [-pi pi], 'periodic');
try
    h = circconv(f, g);
catch ME
    pass(4) = strcmpi(ME.message, 'Domains of f and g must match.');
end

end
