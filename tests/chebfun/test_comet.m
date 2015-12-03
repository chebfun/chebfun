% Test for @chebfun/comet.m.

function pass = test_comet(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Real scalar functions.
fsr1 = chebfun(@sin, [-1 0 1], pref);
fsr2 = chebfun(@cos, [-1 0.5 1], pref);

% Obviously, we can't check if the plots are correct without human intervention,
% so all these tests are meant to do is make sure that we don't crash.

hfig = figure('Visible', 'off');

% Check COMET and COMET3.
pass(1) = doesNotCrash(@() comet(fsr1));
pass(2) = doesNotCrash(@() comet(fsr1, fsr2));

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
