function pass = test_isempty(varargin)

f = funcheb1();
pass(1) = isempty(f);

f = funcheb1(@sin);
pass(2) = ~isempty(f);

f = funcheb1(@(x) [sin(x), cos(x)]);
pass(3) = ~isempty(f);

f = [ funcheb1(@sin), funcheb1(@sin) ];
pass(4) = ~isempty(f);

end
