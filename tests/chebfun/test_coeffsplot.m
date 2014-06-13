% Tests for @chebfun/coeffsplot.m.

function pass = test_coeffsplot(pref)

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

pass(1) = doesNotCrash(@() coeffsplot(f));
pass(2) = doesNotCrash(@() coeffsplot(g));
pass(3) = doesNotCrash(@() coeffsplot(F));
pass(4) = doesNotCrash(@() coeffsplot(G));
pass(5) = doesNotCrash(@() coeffsplot(Q));

% Check plot flags and other options.
pass(6) = doesNotCrash(@() coeffsplot(g, 'noepslevel'));
pass(7) = doesNotCrash(@() coeffsplot(g, 'loglog'));
pass(8) = doesNotCrash(@() coeffsplot(g, '.--'));

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
