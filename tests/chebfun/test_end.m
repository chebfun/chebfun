function pass = test_end(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% test end in infinite dimension:
% Scalar:
f = chebfun(@(x) sin(pi*x));
out = f(end);
pass(1) = isnumeric(out) && length(out) == 1 && abs(out) < epslevel(f);

% Array:
f = chebfun(@(x) [sin(pi*x), cos(pi*x), x]);
out = f(end);
pass(2) = isnumeric(out) && all(size(out) == [1, 3]) && ... 
    norm(out - [0 -1 1], inf) < epslevel(f);

out = f(end,2);
pass(3) = isnumeric(out) && numel(out) == 1 && abs(out + 1) < epslevel(f);

% transpose
g = f.';
out = g(end);
pass(4) = isnumeric(out) && all(size(out) == [3, 1]) && ... 
    norm(out - [0 -1 1].', inf) < epslevel(f);

%% test end in scalar dimension:
f = chebfun(@(x) [sin(pi*x), cos(pi*x), x]);
out = f(0,end);
pass(5) = isnumeric(out) && length(out) == 1&& abs(out) < epslevel(f);

g = f.';
out = f(0,end);
pass(6) = isnumeric(out) && length(out) == 1&& abs(out) < epslevel(f);

%% test use of ':'
f = chebfun(@(x) [sin(pi*x), cos(pi*x), x]);
out = f(end, :);
pass(7) = isnumeric(out) && all(size(out) == [1, 3]) && ... 
    norm(out - [0 -1 1], inf) < epslevel(f);

out = f(:, end);
x = chebfun(@(x) x);
err = normest(out - x);
pass(8) = err < epslevel(f);

% Transpose:
g = f.';
out = g(end, :);
err = normest(out - x.');
pass(9) = err < epslevel(f);

out = g(:, end);
pass(10) = isnumeric(out) && all(size(out) == [3, 1]) && ... 
    norm(out - [0 -1 1].', inf) < epslevel(f);

%% Test on singular function:

% Set a domain:
dom = [-2 7];

pow = -0.5;
op = @(x) (x - dom(2)).^pow.*sin(100*x);
f = chebfun(op, dom, 'exps', [0 pow], 'splitting', 'on');
out = f(end);
pass(11) = ( isnumeric(out) ) && all( size(out) == ones(1,2) );

%% Test for function defined on unbounded domain:

% Set a domain:
dom = [0 Inf];

op = @(x) 0.75+sin(10*x)./exp(x);
f = chebfun(op, dom, 'splitting', 'on');
out = f(end);
pass(12) = ( isnumeric(out) ) && all(size(out) == ones(1, 2) ) && ...
    ( abs(out - 0.75) < 1e1*epslevel(f) );

end

function out = normest(f, dom)

% Generate a few random points to use as test values.
seedRNG(6178);
if ( nargin == 1 )
    x = 2 * rand(100, 1) - 1;
else
    x = sum(dom) * rand(10, 1) - dom(1);
end

out = norm(feval(f, x), inf);

end
