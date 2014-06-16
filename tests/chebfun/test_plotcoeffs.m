% Tests for @chebfun/plotcoeffs.m.

function pass = test_plotcoeffs(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

f = chebfun(@sin, pref);
g = chebfun({@sin, @exp}, [-1 0 1], pref);
F = chebfun(@(x) [sin(x) cos(x) exp(x)], pref);
G = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 0 1], pref);
Q = cheb2quasi(F);

% Obviously, we can't check if the plots are correct without human
% intervention, so all these tests are meant to do is make sure none of the
% plotting functions crash.

hfig = figure('Visible', 'off');

pass(1) = doesNotCrash(@() plotcoeffs(f));
pass(2) = doesNotCrash(@() plotcoeffs(g));
pass(3) = doesNotCrash(@() plotcoeffs(F));
pass(4) = doesNotCrash(@() plotcoeffs(G));
pass(5) = doesNotCrash(@() plotcoeffs(Q));

% Check plot flags and other options.
pass(6) = doesNotCrash(@() plotcoeffs(g, 'noepslevel'));
pass(7) = doesNotCrash(@() plotcoeffs(g, 'loglog'));
pass(8) = doesNotCrash(@() plotcoeffs(g, '.--'));

close(hfig);

end

function pass = doesNotCrash(fn)

try
    fn();
    pass = true;
catch ME;
    pass = false;
    rethrow(ME)
end

end
