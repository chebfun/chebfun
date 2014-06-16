% Test for @chebfun/comet.m and @chebfun/comet3.m.

function pass = test_comet(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Real scalar functions.
fsr1 = chebfun(@sin, [-1 0 1], pref);
fsr2 = chebfun(@cos, [-1 0.5 1], pref);
fsr3 = chebfun(@exp, [-1 -0.5 1], pref);

% Real array-valued functions.
far1 = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 0 1], pref);

% Real quasimatrices.
fqr1 = cheb2quasi(far1);

% Obviously, we can't check if the plots are correct without human
% intervention, so all these tests are meant to do is make sure that we don't
% crash.

hfig = figure('Visible', 'off');

% Check COMET and COMET3.
pass(1) = doesNotCrash(@() comet(fsr1));
pass(2) = doesNotCrash(@() comet(fsr1, fsr2));
pass(3) = doesNotCrash(@() comet3(fsr1, fsr2, fsr3));
pass(4) = doesNotCrash(@() comet3(far1));
pass(5) = doesNotCrash(@() comet3(fqr1));

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
