function pass = test_isempty(pref)

if ( nargin == 0 )
    pref = chebfunpref();
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

dom = [-2 7];
pow = -1.64;
f = chebfun(@(x) (x-dom(1)).^pow, dom, 'exps', [pow 0], 'splitting', 'on');
pass(10) = ~isempty(f);

%% Tests for functions defined on unbounded domain:

% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];

op = @(x) (1-exp(-x.^2))./x;
f = chebfun(op, dom);
pass(11) = ~isempty(f);

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2));
f = chebfun(op, dom, 'exps', [2 2]); 
pass(12) = ~isempty(f);

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
f = chebfun(op, dom);
pass(13) = ~isempty(f);

end
