function pass = test_shortPulses(~)
% This test checks whether short forcing pulses influence the solution returned
% from CHEBOP. See #1512.

% Set up a simple linear IVP chebop
t = chebfun('t',[0 2]);
L = chebop(@(t) diff(t) + t,[0 2], 1, []);

% No pulse reference solution:
uNoPulse = L\0;

%% Check that long pulses work:
longPulse = 20*(t>1).*(t<1.2);
uLongPulse = L\longPulse;
pass(1) = all(uLongPulse.domain == [0 1 1.2 2]) && ...
    norm(uNoPulse - uLongPulse) > .1;

%% Check that short pulses work:
shortPulse = 20*(t>1).*(t<1.05);
uShortPulse = L\shortPulse;
pass(2) = all(uShortPulse.domain == [0 1 1.05 2]) && ... 
    norm(uNoPulse - uShortPulse) > .1;

%% Turn off restarting

% Once we turn off restarting, we expect u = uShortPulse (up to a small error).
% This is discussed in #1512, in particular, @nickhale makes a comment that this
% is a MATLAB fault. If this part of the test starts failing, it might be
% because MATLAB are correctly detecting short pulses.

pref = cheboppref();
pref.ivpRestartSolver = false;

uNoRestart = solveivp(L, shortPulse, pref);
pass(3) = all(uShortPulse.domain == [0 1 1.05 2]) && ... 
    norm(uNoPulse - uNoRestart) < 1e-10;

end