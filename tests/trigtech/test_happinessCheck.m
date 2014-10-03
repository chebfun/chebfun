% Test file for trigtech/happinessCheck.m

function pass = test_happinessCheck(pref)

% Get preferences:
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Set the tolerance:
tol = 50*pref.eps;
    
%%
% Test on a scalar-valued function:
x = testclass.trigpts(33);
omega = 8;
f = @(x) sin(omega*pi*x);
g = testclass.make(f(x));
[ishappy, epslevel, tail] = happinessCheck(g, f, [], pref);
pass(1) = tail == 2*omega+1;
pass(2) = ishappy && epslevel < tol;

%%
% Test on an array-valued function:
omega = 7;
f = @(x) [sin(pi*x) cos(floor(omega/2)*pi*x) (sin(omega*pi*x)+cos(omega*pi*x))];
g = testclass.make(f(x));
[ishappy, epslevel, tail] = happinessCheck(g, f, [], pref);
pass(3) = tail == 2*omega+1;
pass(4) = ishappy && all(epslevel < tol);

%%
k = 4*8; % Multiple of four;
m = k/4;
x = testclass.trigpts(k+1);
f = @(x) sin((k+m+1)*pi*x);

% This should be happy, as aliasing fools the happiness test:
pref.sampleTest = 0;
g = testclass.make(f(x));
[ishappy, epslevel, tail] = happinessCheck(g, f, [], pref);
pass(5) = ( ishappy && tail == 2*m+1);

% This should be unhappy, as sampletest fixes things:
pref.sampleTest = 1;
g = testclass.make(f(x));
[ishappy, epslevel, tail] = happinessCheck(g, f, [], pref);
pass(6) = ~ishappy && tail == 33;

end
