function pass = test_realpow(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% These should fail:

x = chebfun('x', pref);
try
    realpow(x, 1/3);
    pass(1) = false;
catch ME
    pass(1) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:realpow:complexRes');
end

try
    realpow(1i*abs(x), 3);
    pass(2) = false;
catch ME
    pass(2) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:realpow:complexRes');
end

%% These should pass:

realpow(x, 1);
pass(3) = true;

realpow(1i*x, 0);
pass(4) = true;

%% Test array-valued:

X = [x x];
try
    realpow(1i*X, 3)
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:realpow:complexRes');
end

try
    realpow(X,.42);
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:realpow:complexRes');
end

realpow(X, 1);
pass(7) = true;

realpow(1i*X, 0);
pass(8) = true;

end
