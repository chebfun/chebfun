function pass = test_isempty(varargin)

f = funcheb2();
pass(1) = isempty(f);

f = funcheb2(@sin);
pass(2) = ~isempty(f);

f = funcheb2(@(x) [sin(x), cos(x)]);
pass(3) = ~isempty(f);

f = [ funcheb2(@sin), funcheb2(@sin) ];
pass(4) = ~isempty(f);

end
