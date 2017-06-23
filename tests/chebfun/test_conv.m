% Test file for @chebfun/conv.m.

function pass = test_conv(pref)

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end

% Construct a few CHEBFUN objects for the tests.
f = chebfun(@(x) x, [-1 1]);
g = chebfun(@(x) sin(5*x), [2 4]);
h = chebfun(@(x) cos(2*x), [-3 1]);
p = chebfun(@(x) sin(15*x), [-1 1]);
q = chebfun(@(x) exp(cos(3*x)), [-1 1]);

j = 0;

for algo = ["cheb", "leg"]
    
    %% 1. Test the commutativity.
    H1 = conv(f, g, algo);
    H2 = conv(g, f, algo);
    j = j+1;
    pass(j) = norm(H1 - H2) < 100*eps;
    
    %% 2. Test the associativity.
    H1 = conv(conv(f, g, algo), h, algo);
    H2 = conv(f, conv(g, h, algo), algo);
    j = j+1;
    pass(j) = norm(H1 - H2) < 100*eps;
    
    %% 3. Test the distributivity.
    H1 = conv(f, (p + q), algo);
    H2 = conv(f, p, algo) + conv(f, q, algo);
    j = j+1;
    pass(j) = norm(H1 - H2) < 100*eps;
    
    %% 4. Test B-splines.
    % Generates piecewise polynomial cardinal B-splines by a process of
    % convolution, with CONV, then estimates the error compared to the exact
    % solution B, using NORM. For more on generating B-splines using convolution,
    % see e.g. http://en.wikipedia.org/wiki/B-spline
    B = (1/6)*chebfun( {@(x) (2+x).^3, ...
        @(x)1 + 3*(1+x) + 3*(1+x).^2 - 3*(1+x).^3, ...
        @(x)1 + 3*(1-x) + 3*(1-x).^2 - 3*(1-x).^3, ...
        @(x)(2-x).^3}, -2:2 );
    s = chebfun(1, [-.5 .5]);
    f = s;
    for k = 1:3
        f = conv(f, s, algo);
    end
    j = j+1;
    pass(j) = norm(f - B) < 100*eps*vscale(f);
    
    %% 5. Test more splines.
    f = chebfun(1/2);
    g = f;
    for k = 2:4
        g = conv(f, g, algo);
    end
    j = j+1;
    pass(j) = abs(feval(g, 1) - 23/96) < 100*eps;
    
    %% 6. Test the example from the HELP text:
    f = chebfun(1/2); g = f;
    for k = 2:4
        g = conv(f, g, algo);
    end
    g1 = g(.1);
    for k = 1:4
        g = diff(g);
    end
    g2 = g(.1);
    err = abs(g1 - 0.332114583333333) + abs(g2 - 0);
    j = j+1;
    pass(j) = err < 1e-14;
    
    %% 7. Test the example from ATAP
    a = 0.2885554757; b = 0.3549060246;
    g = chebfun(@(x) sin(1./x).*sin(1./sin(1./x)), [a,b], 80, 'chebkind', 2);
    t = 1e-7;
    f = chebfun(@(x) exp(-x.^2/(4*t))/sqrt(4*pi*t),.003*[-1 1]);
    h = conv(f, g, algo);
    err = abs(norm(h,2) - 0.029781437647379);
    tol = 10*max(get(f, 'vscale')*eps, ...
        eps*get(g, 'vscale'));
    j = j+1;
    pass(j) = err < tol;
    
    %% An Example due to Mohsin Javed:
    f = chebfun(@(x) x, [-1.5 0] );
    g = chebfun(@(x) sin(x), [-1 1]);
    h1 = conv(f, g, algo);
    h2 = conv(f, g, 'old');
    tol = 10*max(get(f, 'vscale')*eps, ...
        eps*get(g, 'vscale'));
    j = j+1;
    pass(j) = norm(h1 - h2, inf) < tol;
    
    %% Testing Delta function convolution
    % Delta funciton reproduces the function under convolution:
    f = chebfun({@(x) sin(x), @(x) cos(x), @(x) sin(4*x.^2)}, [-2, -1, 0, 1] );
    x = chebfun('x', [-2, 1] );
    d = dirac(x);
    g = conv(d, f, algo);
    g = restrict(g, [-2, 1]);
    g.pointValues(1) = 2*g.pointValues(1);
    g.pointValues(end) = 2*g.pointValues(end);
    j = j+1;
    pass(j) = norm(f - g, inf) < tol;
    
    % Derivative of delta function differentiates the function:
    x = chebfun('x');
    f = sin(x);
    g = conv(f, diff(dirac(x)), algo);
    g = restrict(g, [-1, 1] );
    g.pointValues(1) = 2*g.pointValues(1);
    g.pointValues(end) = 2*g.pointValues(end);
    j = j+1;
    pass(j) = norm(g - cos(x), inf ) < 10*tol;
    
    % Second order ODE via delta functions and convolutions:
    % g = f'' + f
    f = sin(x);
    g = conv(f, diff(dirac(x), 2) + dirac(x), algo);
    g = restrict(g, [-1, 1] );
    g.pointValues(1) = 2*g.pointValues(1);
    g.pointValues(end) = 2*g.pointValues(end);
    j = j+1;
    pass(j) = norm(g , inf ) < 1e3*tol;
    
    %% Maurice's Cox examples:
    fX1 = chebfun(@(x) exp(x), [0, log(2)]);
    fX2 = chebfun(@(x) exp(x), [log(2), log(3)]);
    fX3 = chebfun(@(x) exp(x), [log(3), log(4)]);
    g1 = conv(fX1, fX3, algo);
    g2 = conv(g1, fX2, algo);
    
    g3 = conv(fX1, fX2, algo);
    g4 = conv(g3, fX3, algo);
    j = j+1;
    pass(j) = normest( g2 - g4 ) < 1e1*tol;
    
    g5 = conv(fX2, fX3, algo);
    g6 = conv(g5, fX1, algo);
    j = j+1;
    pass(j) = normest( g2 - g6 ) < 1e1*tol;
    
    %% test 'same' option
    f = chebfun(@(x) exp(-x.^2), [-10 10]);
    g = chebfun(@(x) exp(-x.^2), [-20 20]);
    h = conv(f, f, 'same', algo);
    j = j+1;
    pass(j) = norm(h.domain([1, end]) - [-10 10], inf) < eps*10;
    
    h = conv(f, g, 'same', algo);
    j = j+1;
    pass(j) = norm(h.domain([1, end]) - [-10 10], inf) < eps*10;
    
    %% test quasimatrix option
    f = chebfun(@sin);
    g = chebfun(@cos);
    fg = [f g];
    gg = [g g];
    h1 = conv(fg, gg, algo);
    h2 = [conv(f,g, algo), conv(g,g, algo)];
    j = j+1;
    pass(j) = norm(h1 - h2) < eps*10;
    
    ffg = conv([f f],g, algo);
    fgg = conv(f,[g g], algo);
    j = j+1;
    pass(j) = norm(ffg-fgg) < 10*eps;
    
    cheb.x;
    phi = @(t) chebfun(@(x) exp(-x^2/(4*t))/sqrt(4*pi*t));
    f = 1 - 0.2*x - abs(x-0.2);
    fsmooth = conv(f,phi(1e-4),'same', algo);
    j = j+1;
    pass(j) = norm(f-fsmooth) < .01;
    
end

end

