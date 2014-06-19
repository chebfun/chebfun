% Test file for chebfun trigonometric and related functions.

function pass = test_trig(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

pref.splitting = 1;

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% List of trigonometric functions to test.
trigFunctions = {@acos, @acosd, @acosh, @acot, @acotd, @acoth, @acsc, ...
    @acscd, @acsch, @asec, @asecd, @asech, @asin, @asind, @asinh, @atan, ...
    @atand, @atanh, @cos, @cosd, @cosh, @cot, @cotd, @coth, @csc, @cscd, ...
    @csch, @sec, @secd, @sech, @sin, @sind, @sinh, @tan, @tand, @tanh, @mysinc};

% [TODO]: Add tests for ATAN2() and ATAN2D().

% Function with which we will be composing.  (The reason for the shift and
% scaling, etc. is to prevent any problems with the functions under test from
% evaluating to Inf.)
base_op = @(x) sign(x - 0.1).*abs(x + 0.2).*sin(3*x)*(pi/16) + pi/8;
f = chebfun(base_op, [-1 -0.2 0.1 1], pref);

% Do the tests.
for k = 1:1:numel(trigFunctions)
    trig_op = trigFunctions{k};
    g_exact = @(x) trig_op(base_op(x));
    g = trig_op(f, pref);
    err = feval(g, xr) - g_exact(xr);
    pass(k) = norm(err, inf) < 50*g.vscale.*g.epslevel;
end

end

function y = mysinc(x, varargin)
% Chebfun SINC differs in definition from MATLAB's SINC, so define a special
% function for it here.
    y = sin(x)./x;
end

