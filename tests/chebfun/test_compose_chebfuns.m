% Test file for @chebfun/compose.m (composition of chebfuns).

function pass = test_compose_chebfuns(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Create preference structure with splitting enabled.
pref_split = pref;
pref_split.splitting = 1;

% [TODO]:  Once composeChebfuns() can handle domain checks, test that too.

% Two smooth chebfuns.
pass(1) = test_one_compose_chebfuns(@(x) cos(2*(x + 0.2)), ...
    @(x) sin(x - 0.1), [-1 1], pref);

% Two smooth chebfuns, non-default domain.
pass(2) = test_one_compose_chebfuns(@(x) cos(2*(x + 0.2)), ...
    @(x) sin(x - 0.1), [-2 7], pref);

% Smooth chebfun of a non-smooth chebfun.
pass(3) = test_one_compose_chebfuns(@(x) abs(cos(2*(x + 0.2))), ...
    @(x) sin(x - 0.1), [-1 1], pref_split);

% Non-smooth chebfun of a smooth chebfun.
pass(4) = test_one_compose_chebfuns(@(x) sin(x - 0.1), ...
    @(x) abs(cos(2*(x + 0.2))), [-1 1], pref_split);

% Non-smooth chebfun of a non-smooth chebfun.
pass(5) = test_one_compose_chebfuns(@(x) abs(sin(2*(x - 0.1))), ...
    @(x) abs(cos(2*(x + 0.2))), [-1 1], pref_split);

% Single-column chebfun of an array-valued chebfun.
pass(6) = test_one_compose_chebfuns(@(x) [sin(x - 0.1) cos(x - 0.2)], ...
    @(x) exp(x), [-1 1], pref);

% Array-valued chebfun of a single-column chebfun.
pass(7) = test_one_compose_chebfuns(@(x) exp(x)/exp(1), ...
    @(x) [sin(x - 0.1) cos(x - 0.2)], [-1 1], pref);

% Single-column chebfun of an array-valued chebfun.
f_exact = @(x) [sin(x - 0.1) cos(x - 0.2)];
g_exact = @(x) exp(x);
fq = quasimatrix(f_exact, [-1 1], pref);
g = chebfun(g_exact, [-1 1], pref);
pass(8) = test_one_compose_chebfuns_quasi(f_exact, g_exact, [-1 1], pref, fq, g);

% % Array-valued chebfun of a single-column chebfun.
f_exact = @(x) exp(x)/exp(1);
g_exact = @(x) [sin(x - 0.1), cos(x - 0.2)];
f = chebfun(f_exact, [-1 1], pref);
gq = quasimatrix(g_exact, [-1 1], pref);
pass(9) = test_one_compose_chebfuns_quasi(f_exact, g_exact, [-1 1], pref, f, gq);

% Can't compose two array-valued chebfuns.
try
    test_one_compose_chebfuns(@(x) [sin(x) cos(x)], @(x) [exp(x) -sin(x)], ...
        [-1 1], pref);
    pass(10) = false;
catch ME
    pass(10) = true;
end

%% Compose a function defined on unbounded domain with a function defined on
% bounded domain, i.e. G(F):

% Set the domain:
dom = [0 Inf];
domCheck = [0 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

opf = @(x) exp(-x);
opg = @(x) cos(x);
oph = @(x) cos(exp(-x));
f = chebfun(opf, dom);
g = chebfun(opg, [-1 1]);
h = compose(f, g);
hVals = feval(h, x);
hExact = oph(x);
err = hVals - hExact;
pass(11) = norm(err, inf) < 1e1*eps*get(h,'vscale');


%% Test compose with a discontinuous breakpoint. See #1074.
ep = 0.25;
x = chebfun('x');
F = (abs(x) < ep)/(2*ep);
pass(12) = get(F(x), 'ishappy');


%% Test composition g(f) of a periodic chebfun f and a chebfun g.  See #1963.
f = chebfun(@(x) sin(pi*x), 'trig');
cheb.x
g = x.^2;
h = g(f);
pass(13) = isPeriodicTech(h);

G = [exp(x), abs(x)];
H = G(f);
pass(14) = isPeriodicTech(H(:,1));
pass(15) = ~isPeriodicTech(H(:,2));

%% Test composition with a CHEBFUN2:
Fr = chebfun(@(t) [ cos(t), sin(t) ], [ -pi, pi ]); % Two real columns.
Fc = chebfun(@(t) exp(1i*t), [ -pi, pi ]);          % One complex column.
g = chebfun2(@(x,y) x.^2 + y.^2);
h_true = chebfun(@(t) 1 + 0*t, [ -pi, pi ]);
hr = compose(Fr, g);
hc = compose(Fc, g);
pass(16) = ( norm(hr - h_true) < 100 * eps );
pass(17) = ( norm(hc - h_true) < 100 * eps );

%% Test composition of a periodic CHEBFUN with a CHEBFUN2:
F = chebfun(@(t) [ cos(t), sin(t) ], [ -pi, pi ], 'trig');
g = chebfun2(@(x,y) x + y);
h_true = chebfun(@(t) cos(t) + sin(t), [ -pi, pi ], 'trig');
h = compose(F, g);
pass(18) = ( norm(h - h_true) < 100 * eps );
pass(19) = isPeriodicTech(h);

%% Test composition [periodic, nonperiodic] with a CHEBFUN2:
F = [ chebfun(@(t) cos(t), [ -pi, pi ], 'trig'), ...
    chebfun(@(t) sin(t), [ -pi, pi ]) ];
g = chebfun2(@(x,y) x + y);
h = compose(F, g);
pass(20) = ~isPeriodicTech(h);

%% Test composition with a CHEBFUN2V:
Fr = chebfun(@(t) [ cos(t), sin(t) ], [ -pi, pi ]); % Two real columns.
Fc = chebfun(@(t) exp(1i*t), [ -pi, pi ]);          % One complex column.
G = chebfun2v(@(x,y) x.^2 + y.^2, @(x,y) x);
H_true = chebfun(@(t) [ 1 + 0*t, cos(t) ], [ -pi, pi ]);
Hr = compose(Fr, G);
Hc = compose(Fc, G);
pass(21) = ( norm(Hr - H_true) < 100 * eps );
pass(22) = ( norm(Hc - H_true) < 100 * eps );

%% Test composition of a periodic CHEBFUN with a CHEBFUN2V
F = chebfun(@(t) [ cos(t), sin(t) ], [ -pi, pi ], 'trig');
g = chebfun2v(@(x,y) x, @(x,y) y);
h_true = chebfun(@(t) [ cos(t) , sin(t) ], [ -pi, pi ], 'trig');
h = compose(F, g);
pass(23) = ( norm(h - h_true) < 100 * eps );
pass(24) = isPeriodicTech(h);

%% Test composition with a CHEBFUN3:
F = chebfun(@(t) [ cos(t), sin(t), t ], [ -pi, pi ]);
g = chebfun3(@(x,y,z) x.^2 + y.^2 + z.^2);
h_true = chebfun(@(t) 1 + t.^2, [ -pi, pi ]);
h = compose(F, g);
pass(25) = ( norm(h - h_true) < 100 * eps );

F = chebfun(@(t) [ cos(t), sin(t), 1i*t ], [ -pi, pi ]);
g = chebfun3(@(x,y,z) x.^2 + y.^2 + z.^2);
try
    h = compose(F, g);
    pass(26) = 0;
catch EM
   pass(26) = strcmp(EM.identifier, 'CHEBFUN:CHEBFUN:compose:complex3');
end

%% Test composition with a CHEBFUN3V:
F = chebfun(@(t) [ cos(t), sin(t), t ], [ -pi, pi ]);
G = chebfun3v(@(x,y,z) x.^2 + y.^2, @(x,y,z) z);
H_true = chebfun(@(t) [ 1 + 0*t, t ], [ -pi, pi ]);
H = compose(F, G);
pass(27) = ( norm(H - H_true) < 100 * eps );

%% Test composition of a periodic CHEBFUN with a CHEBFUN3 and CHEBFUN3V:
F = chebfun(@(t) [ cos(t), sin(t), cos(2*t) ], [ -pi, pi ], 'trig');
g = chebfun3(@(x,y,z) x.^2 + y.^2 + z.^2);
h_true = chebfun(@(t) 1 + cos(2*t).^2, [ -pi, pi ], 'trig');
h = compose(F, g);
pass(28) = ( norm(h - h_true) < 10 * eps );
pass(29) = isPeriodicTech(h);

G = chebfun3v(@(x,y,z) x.^2 + y.^2, @(x,y,z) z);
H = G(F);
pass(30) = isPeriodicTech(H);


end

% Test composition of two chebfuns.
function pass = test_one_compose_chebfuns(f_exact, g_exact, dom, pref)
    % Random points to use as test values.
    seedRNG(7681);
    xr = 2 * rand(100, 1) - 1;

    f = chebfun(f_exact, dom, pref);
    g = chebfun(g_exact, dom, pref);
    h = compose(f, g, [], pref);
    h_exact = @(x) g_exact(f_exact(x));
    x = ((dom(end) - dom(1))/2)*xr + dom(1) + (dom(end) - dom(1))/2;
    err = norm(feval(h, x) - h_exact(x), inf);
    pass = (err < 20*vscale(h)*eps) && ...
        isequal(feval(g, feval(f, dom)), feval(h, dom));
end

% Test composition of two chebfuns or quasimatrices.
function pass = test_one_compose_chebfuns_quasi(f_exact, g_exact, dom, pref, f, g)
    % Random points to use as test values.
    seedRNG(7681);
    xr = 2 * rand(100, 1) - 1;
    h = compose(f, g, [], pref);
    h_exact = @(x) g_exact(f_exact(x));
    x = ((dom(end) - dom(1))/2)*xr + dom(1) + (dom(end) - dom(1))/2;
    err = norm(feval(h, x) - h_exact(x), inf);
    pass = (err < 20*vscale(h)*eps) && ...
        isequal(feval(g, feval(f, dom)), feval(h, dom));
end
