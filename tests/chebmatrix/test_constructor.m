function pass = test_constructor(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Test construction from a matrix of values;
seedRNG(42)
r = rand(2);
A = chebmatrix(r);
B = chebmatrix(num2cell(r));
pass(1) = norm(A - B) < 1e-14;

% This shouldn't work:
try
    C = chebmatrix({r});
catch ME
    pass(2) = strcmp(ME.identifier, 'CHEBFUN:CHEBMATRIX:input:nonscalarcell');
end

x = chebfun('x');
f = [sin(x) cos(x)];

% This shouldn't work:
try
    D = chebmatrix({f});
catch ME
    pass(3) = strcmp(ME.identifier, 'CHEBFUN:CHEBMATRIX:input:arraychebfun');
end

% Chebmatrix with num2cell the inputs:
E = chebmatrix(f);
pass(4) = numel(E.blocks) == 2;

end

