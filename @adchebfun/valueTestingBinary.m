function err = valueTestingBinary(func)

% TODO: Document.

%% Initialise

% Seed random generator to ensure same values.
seedRNG(6179);

% Length of test functions:
N = 8;

% Generate two arbitrary CHEBFUN objects to evaluate the function at:
u1 = chebfun(0.1*rand(N, 1) + .5);
u2 = chebfun(0.1*rand(N,1) + .5);

% Construct corresponding ADCHEBFUN objects:
v1 = adchebfun(u1);
v2 = adchebfun(u2);

% Also create two arbitrary CHEBFUNS to test with:
w1 = chebfun(0.1*rand(N, 1) + .5);
w2 = chebfun(0.1*rand(N,1) + .5);

% And finally two scalars to test with:
s1 = rand();
s2 = rand();

% Initialise error vector:
err = zeros(1,5);

%% Create various combinations

% ADCHEBFUN and ADCHEBFUN:
v1v2 = func(v1, v2);

% ADCHEBFUN and CHEBFUN:
v1w2 = func(v1, w2);
w1v2 = func(w1, v2);

% ADCHEBFUN and SCALAR:
v1s2 = func(v1, s2);
s1v2 = func(s1, v2);

% Do CHEBFUN counterparts for allowing comparing values:
u1u2 = func(u1, u2);
u1w2 = func(u1, w2);
w1u2 = func(w1, u2);
u1s2 = func(u1, s2);
s1u2 = func(s1, u2);

%% Compare values
err(1) = norm(v1v2.func - u1u2);
err(2) = norm(v1w2.func - u1w2);
err(3) = norm(w1v2.func - w1u2);
err(4) = norm(v1s2.func - u1s2);
err(5) = norm(s1v2.func - s1u2);

% Check correctness of derivatives
%     v1v2der = v1v2.jacobian;
%     v1w2der = v1w2.jacobian;
%     w1v2der = w1v2.jacobian;
%     v1s2der = v1s2.jacobian;
%     s1v2der = s1v2.jacobian;

%     % Compare derivatives with expected values. For third and fifth operation,
%     % need to call func with 0 as the first argument to get the correct sign.
%     pass(funcCounter + 2, 1) = ...
%         ( norm(v1v2der*w1 - (func(v1.jacobian, v2.jacobian))*w1 ) == 0 );
%     pass(funcCounter + 2, 2) = ...
%         ( norm(v1w2der*w2 - v1.jacobian*w2 ) == 0 );
%     pass(funcCounter + 2, 3) = ...
%         ( norm(w1v2der*w1 - func(0,v2.jacobian)*w1 ) == 0 );
%     pass(funcCounter + 2, 4) = ...
%         ( norm(v1s2der*w2 - v1.jacobian*w2 ) == 0 );
%     pass(funcCounter + 2, 5) = ...
%         ( norm(s1v2der*w1 - func(0,v2.jacobian)*w1 ) == 0 );


end