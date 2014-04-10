% Test for linearity detection of various complicated expresssions

function pass = test_linearityDetection

%% Initialise
% Seed random generator to ensure same values.
seedRNG(6179);

% Length of test functions:
N = 8;

% Generate two arbitrary CHEBFUN objects to evaluate the function at:
u1 = chebfun(0.1*rand(N, 1) + .5);
u2 = chebfun(0.1*rand(N,1) + .5);

% Construct corresponding ADCHEBFUN objects. Here, we seed the ADCHEBFUN objects
% to ensure that the dimensions of their derivatives are correct:
v1 = seed(adchebfun(u1), 1, 2);
v2 = seed(adchebfun(u2), 2, 2);

% Also create two arbitrary CHEBFUNS to test with:
w1 = chebfun(0.1*rand(N, 1) + .5);
w2 = chebfun(0.1*rand(N,1) + .5);

% And finally two scalars to test with:
s1 = rand();
s2 = rand();

%% Check linearity of various expressions
pass = [];

f = sin(v1) + w1;
pass(length(pass) + 1) = all( f.linearity == [0 1]);

f = sin(v1) + v1;
pass(length(pass) + 1) = all( f.linearity == [0 1]);

f = cos(v1) + exp(w2);
pass(length(pass) + 1) = all( f.linearity == [0 1]);

f = log(v2 + 2) + exp(s1);
pass(length(pass) + 1) = all( f.linearity == [1 0]);

f = sin(v1) - exp(v2);
pass(length(pass) + 1) = all( f.linearity == [0 0]);

f = sin(v1) - diff(v2);
pass(length(pass) + 1) = all( f.linearity == [0 1]);

f = sin(v1) - diff(v2);
pass(length(pass) + 1) = all( f.linearity == [0 1]);

f = tan(v2).*(w2);
pass(length(pass) + 1) = all( f.linearity == [1 0]);

f = w2.*log2(v1+2)/4;
pass(length(pass) + 1) = all( f.linearity == [0 1]);

%% Three ADCHEBFUN objects involved
u3 = chebfun(0.1*rand(N,1) + .5);

v1 = seed(adchebfun(u1), 1, 3);
v2 = seed(adchebfun(u2), 2, 3);
v3 = seed(adchebfun(u3), 3, 3);

f = v1.*v3 + v2;
pass(length(pass) + 1) = all( f.linearity == [0 1 0]);

f = v1./(v3+2) + v2;
pass(length(pass) + 1) = all( f.linearity == [0 1 0]);

f = v1 + v2.^2 + csc(v3);
pass(length(pass) + 1) = all( f.linearity == [1 0 0]);

end