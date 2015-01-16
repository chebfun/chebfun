% Test file for trigtech/simplify.m

function pass = test_simplify(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

testclass = trigtech();

% Tolerance for passing to simplify:
simptol = 1e-6;

%%
% Test pathological inputs.

% Empty TRIGTECH objects should be left alone.
f = testclass.make();
g = simplify(f);
pass(1) = isequal(f, g);

% Unhappy TRIGTECH objects should be left alone.
f = testclass.make(@(x) sin(x));
g = simplify(f);
pass(2) = ~f.ishappy && isequal(f, g);

%%
% Test for a scalar-valued function:

f = testclass.make(@(x) exp(sin(2*pi*x)) + exp(cos(3*pi*x)));
g = simplify(f, simptol);
pass(3) = abs(g.coeffs(end)) ~= 0;
pass(4) = length(g) < length(f);
pass(5) = norm(feval(f, x) - feval(g, x), inf) < 10*g.epslevel*g.vscale;

%%
% Lengths of simplifications should be invariant under scaling:

f1 = 1e-8*f;
g1 = simplify(f1, simptol);
pass(6) = length(g1) == length(g);

f2 = 1e8*f;
g2 = simplify(f2, simptol);
pass(7) = length(g2) == length(g);

%%
% Test for an array-valued function:

f = testclass.make(@(x) [exp(sin(2*pi*x)) exp(cos(3*pi*x)) 3./(4-cos(pi*x))]);
g = simplify(f, simptol);
pass(8) = any(abs(g.coeffs(1, :)) ~= 0);
pass(9) = length(g) < length(f);
pass(10) = all(norm(feval(f, x) - feval(g, x), inf) < ...
    10*max(g.epslevel.*g.vscale));

%%
% Try a contrived example which will leave a zero TRIGTECH:

f = testclass.make(@(x) sin(100*pi*(x + 0.1)));
g = simplify(f, 1e20);
pass(11) = iszero(g);

%%
% Check that a long identically-zero TRIGTECH simplifies correctly:
f = testclass.make(@(x) 0*x, [], struct('fixedLength', 8));
g = simplify(f);
pass(12) = iszero(g) && (length(g) == 1);

end
