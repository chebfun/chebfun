% Tests for @chebfun/chebpolyplot.m.

function pass = test_chebpolyplot(pref)

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

pass(1) = doesNotCrash(@() chebpolyplot(f));
pass(2) = doesNotCrash(@() chebpolyplot(g));
pass(3) = doesNotCrash(@() chebpolyplot(F));
pass(4) = doesNotCrash(@() chebpolyplot(G));
pass(5) = doesNotCrash(@() chebpolyplot(Q));

% Check plot flags and other options.
pass(6) = doesNotCrash(@() chebpolyplot(g, 'noepslevel'));
pass(7) = doesNotCrash(@() chebpolyplot(g, 'loglog'));
pass(8) = doesNotCrash(@() chebpolyplot(g, '.--'));

close(hfig);

end

function pass = doesNotCrash(fn)

try
    fn();
    pass = true;
catch ME;
    pass = false;
end

end
