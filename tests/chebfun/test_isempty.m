function pass = test_isempty(pref)

if ( nargin == 0 )
    pref = chebpref();
end

f = chebfun();
pass(1) = isempty(f);

f = chebfun([]);
pass(2) = isempty(f);

f = chebfun(@sin, 0);
pass(3) = isempty(f);

f = chebfun([], [-2, 2], pref);
pass(4) = isempty(f);

f = chebfun([], [pi, 42], pref);
pass(5) = isempty(f);

f = chebfun(@sin);
pass(6) = ~isempty(f);

f = chebfun(@sin, [-2 2], pref);
pass(7) = ~isempty(f);

f = chebfun(@sin, [-2 2], 1);
pass(8) = ~isempty(f);

x = chebfun('x');
f = [x abs(x)];
pass(9) = ~isempty(f);

end

