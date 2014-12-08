function pass = test_trigcasting(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfunpref();
end

% Testing of various casting between TRIGFUNS and CHEBFUNS.
% NOTE: The basic operations such as CHEBFUN/PLUS, CHEBFUN/TIMES, 
% CHEBFUN/RESTRICT and CHEBFUN/HORZCAT are also tested in their corresponding 
% tests.

% We are testing here: abs, cumsum, power, round, ceil, floor, quasimatrices, 
% array-valued.
f = chebfun(@(x) cos(pi*x).^2, [1, 3], 'trig' );
pass(1) = abs(sum(f{2,3}) - 1/2) < vscale(f)*epslevel(f)*100;
pass(2) = abs(f(20)-1) < vscale(f)*epslevel(f)*100;

%% Test abs and quasimatrices
f = chebfun(@(x) sin(10*x), [0, 2*pi], 'trig' );
g = chebfun(@(x) x, [0, 2*pi]);
H = [f, g];
tech = pref.tech;
pass(3) = isequal(get(H(:, 1).funs{1}, 'tech'), @trigtech );
pass(4) = isequal(get(H(:, 2).funs{1}, 'tech'), tech );
AH = abs(H);
pass(5) = isequal(get(AH(:, 1).funs{1}, 'tech'), tech );
pass(6) = isequal(get(AH(:, 2).funs{1}, 'tech'), tech );

pass(7) = norm(AH(:,1) - abs(f), inf) < 100*epslevel(f)*vscale(f);
pass(8) = norm(AH(:,2) - abs(g), inf) < 100*epslevel(f)*vscale(f);

%% Test conversion to array-valued
G = quasi2cheb(H);
pass(9)  = isequal(get(G(:, 1).funs{1}, 'tech'), tech);
pass(10) = isequal(get(G(:, 2).funs{1}, 'tech'), tech);
pass(11) = length(G(:,1)) == length(G(:,2));

%% Test for times
f = chebfun(@(x) exp(-x.^2), [-10, 10], 'trig' );
g = chebfun(@(x) exp(-x.^2), [-10, 10]);
x = chebfun(@(x) x, [-10, 10] );
xx = -3:7;
fx = f.*x;
gx = g.*x;
pass(12) = norm(fx(xx)-gx(xx), inf) < 100*epslevel(fx)*vscale(fx);
pass(13) = isequal(get(fx.funs{1}, 'tech'), tech);

%% Test for round, floor, ceil
f = chebfun(@(x) exp(sin(x)), [0 2*pi], 'trig');
g = chebfun(@(x) exp(sin(x)), [0 2*pi]);
pass(14) = norm(round(f) - round(g), inf) < 100*epslevel(f)*vscale(f);
pass(15) = norm(floor(f) - floor(g), inf) < 100*epslevel(f)*vscale(f);
pass(16) = norm(ceil(f) - ceil(g), inf) < 100*epslevel(f)*vscale(f);

%% Test for cumsum
f = chebfun(@(x) sin(2*pi*x), 'trig');
g = cumsum(f);
h = cumsum(f+1e-5);
pass(17) = isequal(get(g.funs{1}, 'tech'), @trigtech);
pass(18) = isequal(get(h.funs{1}, 'tech'), tech);

%% Test diff of a quasimatrix
f = chebfun(@(x) sin(2*pi*x), 'trig');
g = chebfun(@(x) sin(2*pi*x));
H = [f, g];
pass(19) = norm(diff(H, 1, 2)) < 100*vscale(f)*epslevel(f);

%% Test for power
f = chebfun(@(x) 2+sin(pi*x), 'trig');
h = chebfun(@(x) 2+sin(pi*x));
g = chebfun(@(x) exp(x));
pass(20) = norm(f.^g - h.^g, inf) < 100*vscale(f)*epslevel(f); 

%% Test for breakpoints
f = chebfun(@(x) sin(2*pi*x), 'trig');
g = chebfun(f, [-1, 0, 1]);
pass(21) = norm(domain(g) - [-1 0 1], inf) < 10*eps;

end
