function pass = test_isHappy(pref)

% Get preferences:
if ( nargin < 1 )
    pref = funcheb2.pref;
end
% Set the tolerance:
tol = 10*pref.funcheb2.eps;

%%
% Test on a scalar-valued function:
x = funcheb2.chebpts(33);
f = @(x) sin(x);
[ishappy, epslevel, tail] = funcheb2.happinessCheck(f, f(x), [], [], pref);
pass(1) = tail == 14;
pass(2) = ishappy && epslevel < tol;

%%
% Test on a vector-valued function:
f = @(x) [sin(x) cos(x) exp(x)];
[ishappy, epslevel, tail] = funcheb2.happinessCheck(f, f(x), [], [], pref);
pass(3) = tail == 15;
pass(4) = ishappy && epslevel < tol;

%%
n = 32;
m = n/2;
x = funcheb2.chebpts(n+1);
f = @(x) cos((2*n+m)*acos(x));

% This should be happy, as aliasing fools the happiness test:
pref.funcheb2.sampletest = 0;
[ishappy, epslevel, tail] = funcheb2.happinessCheck(f, f(x), [], [], pref);
pass(5) = ( ishappy && tail == 17);

% This should be unhappy, as sampletest fixes things:
pref.funcheb2.sampletest = 1;
[ishappy, epslevel, tail] = funcheb2.happinessCheck(f, f(x), [], [], pref);
pass(6) = ~ishappy && tail == 33;

end