function pass = test_happinessCheck(pref)

% Get preferences:
if ( nargin < 1 )
    pref = funcheb1.pref;
end
% Set the tolerance:
tol = 10*pref.funcheb1.eps;

%%
% Test on a scalar-valued function:
x = funcheb1.chebpts(33);
f = @(x) sin(x);
g = funcheb1(f(x));
[ishappy, epslevel, tail] = happinessCheck(g, f, pref);
pass(1) = tail == 14;
pass(2) = ishappy && epslevel < tol;

%%
% Test on a vector-valued function:
f = @(x) [sin(x) cos(x) exp(x)];
g = funcheb1(f(x));
[ishappy, epslevel, tail] = happinessCheck(g, f, pref);
pass(3) = tail == 15;
pass(4) = ishappy && epslevel < tol;

%%
n = 32;
m = n/2;
x = funcheb1.chebpts(n+1);
f = @(x) cos((2*n+m)*acos(x));

% This should be happy, as aliasing fools the happiness test:
pref.funcheb1.sampletest = 0;
g = funcheb1(f(x));
[ishappy, epslevel, tail] = happinessCheck(g, f, pref);
pass(5) = ( ishappy && tail == 15);

% This should be unhappy, as sampletest fixes things:
pref.funcheb1.sampletest = 1;
g = funcheb1(f(x));
[ishappy, epslevel, tail] = happinessCheck(g, f, pref);
pass(6) = ~ishappy && tail == 33;

end