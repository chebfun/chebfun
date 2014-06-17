% Test file for @chebfun/realsqrt.m.

function pass = test_realsqrt(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% These should fail:

x = chebfun('x', pref);
try
    realsqrt(x);
    pass(1) = false;
catch ME
    pass(1) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:realsqrt:complexRes');
end

try
    realsqrt(1i*abs(x));
    pass(2) = false;
catch ME
    pass(2) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:realsqrt:complexRes');
end

%% These should pass:

realsqrt(x.^2);
pass(3) = true;

realsqrt(abs(x));
pass(4) = true;

%% Test array-valued:

X = [x x];
try
    realsqrt(1i*X);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:realsqrt:complexRes');
end

try
    realsqrt(X);
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:realsqrt:complexRes');
end

realsqrt(X.^2);
pass(7) = true;

realsqrt(abs(X));
pass(8) = true;

end
