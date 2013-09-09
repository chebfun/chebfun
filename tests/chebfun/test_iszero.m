function pass = test_iszero(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

% Test scalars:
f = chebfun(0, pref);
pass(1) = iszero(f);

f = chebfun([], pref);
pass(2) = iszero(f);

f = chebfun(2, pref);
pass(3) = ~iszero(f);

% Test piecewise domains:
f = chebfun(0, [-1, 0, 1], pref);
pass(4) = iszero(f);

f = chebfun([], [-1, 0, 1], pref);
pass(5) = iszero(f);

f = chebfun({0, 1}, [-1, 0, 1], pref);
pass(6) = ~iszero(f);

f = chebfun({1, 0}, [-1, 0, 1], pref);
pass(7) = ~iszero(f);

% Test arrays:
f = chebfun([0 0], pref);
pass(8) = all(iszero(f));

f = chebfun([0 1], pref);
pass(9) = all(iszero(f) == [1 0]);

% [TODO]: Add these once SUBSREF is implemented.
% f = chebfun(0, pref);
% f(0) = 1;
% pass(10) = iszero(f);
% 
% f = chebfun([0 0], pref);
% f(0) = [0, 1];
% pass(11) = all(iszero(f) == [1, 0]);

end